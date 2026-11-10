"""Sensitivity analyses of the Mendelian randomisation estimates (added 2026-11-10, fourth round).

For every selected estimate (gene x study of mr_pqtl.parquet and mr_colon.parquet), from its harmonised instruments:
- heterogeneity: Cochran's Q of the per-instrument Wald ratios (>= 2 instruments);
- MR-Egger (>= 3 instruments): weighted regression of the outcome effect on the exposure effect with an intercept,
  instruments oriented so the exposure effect is positive; the intercept tests directional pleiotropy and the slope
  is the pleiotropy-robust estimate;
- weighted median (>= 3 instruments): median of the Wald ratios weighted by inverse variance, standard error by
  1,000 parametric bootstrap draws (seed 7);
- Steiger direction (all): the variance explained in the exposure must exceed the variance explained in UC, with
  r2 approximated as z2 / (z2 + N) per instrument and summed (N exposure from the Open Targets study index; N outcome
  45,975 for de Lange 2017). A failing direction suggests the variant acts on UC first (reverse causation).
Outputs: data/processed/mr_sensitivity.parquet, mr_sensitivity_summary.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .paths import PROCESSED, RAW

N_OUTCOME = 12366 + 33609
RNG_SEED, N_BOOT = 7, 1000


def _wmedian(b, w):
    o = np.argsort(b)
    b, w = b[o], w[o] / w.sum()
    cum = np.cumsum(w) - 0.5 * w
    return float(np.interp(0.5, cum, b))


def per_estimate(d: pd.DataFrame, n_exp: float, rng) -> dict:
    bx, by, sy = d.beta_exposure.to_numpy(float), d.beta_outcome.to_numpy(float), d.se_outcome.to_numpy(float)
    sx = d.se_exposure.to_numpy(float)
    ratio, se_r = by / bx, sy / np.abs(bx)
    w = 1 / se_r ** 2
    ivw = float((w * ratio).sum() / w.sum())
    out = {"n": len(d)}
    if len(d) >= 2:
        q = float((w * (ratio - ivw) ** 2).sum())
        out.update(cochran_q=q, het_p=float(stats.chi2.sf(q, len(d) - 1)))
    if len(d) >= 3:
        sgn = np.sign(bx)
        X = np.column_stack([np.ones(len(d)), np.abs(bx)])
        W = np.diag(1 / sy ** 2)
        Y = by * sgn
        cov = np.linalg.pinv(X.T @ W @ X)
        coef = cov @ X.T @ W @ Y
        resid = Y - X @ coef
        phi = max(float(resid @ W @ resid) / (len(d) - 2), 1.0)
        se = np.sqrt(np.diag(cov) * phi)
        out.update(egger_slope=float(coef[1]), egger_slope_p=float(2 * stats.t.sf(abs(coef[1] / se[1]), len(d) - 2)),
                   egger_intercept=float(coef[0]), egger_intercept_p=float(2 * stats.t.sf(abs(coef[0] / se[0]), len(d) - 2)))
        wm = _wmedian(ratio, w)
        boots = [_wmedian(rng.normal(by, sy) / rng.normal(bx, sx), w) for _ in range(N_BOOT)]
        wse = float(np.std(boots))
        out.update(wmedian=wm, wmedian_se=wse, wmedian_p=float(2 * stats.norm.sf(abs(wm / max(wse, 1e-12)))))
    zx, zy = bx / sx, by / sy
    r2x = float((zx ** 2 / (zx ** 2 + n_exp)).sum())
    r2y = float((zy ** 2 / (zy ** 2 + N_OUTCOME)).sum())
    out.update(r2_exposure=r2x, r2_outcome=r2y, steiger_correct_direction=bool(r2x > r2y))
    return out


def run() -> dict:
    st = pd.read_parquet(RAW / "extra" / "ot_study.parquet", columns=["studyId", "nSamples"]).set_index("studyId")["nSamples"]
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for est_f, ins_f, layer_col in (("mr_pqtl.parquet", "mr_instruments.parquet", None),
                                    ("mr_colon.parquet", "mr_colon_instruments.parquet", "layer")):
        r = pd.read_parquet(PROCESSED / est_f)
        h = pd.read_parquet(PROCESSED / ins_f)
        hg = {k: d for k, d in h.groupby(["gene", "studyId"])}
        for x in r.itertuples():
            d = hg.get((x.gene, x.studyId))
            if d is None:
                continue
            n_exp = float(st.get(x.studyId, np.nan)) if pd.notna(st.get(x.studyId, np.nan)) else 500.0
            rows.append({"gene": x.gene, "studyId": x.studyId,
                         "layer": getattr(x, layer_col) if layer_col else "plasma_protein",
                         "beta": x.beta, "q": x.q, **per_estimate(d, n_exp, rng)})
    out = pd.DataFrame(rows)
    out.to_parquet(PROCESSED / "mr_sensitivity.parquet", index=False)
    summ = {}
    for layer, d in out.groupby("layer"):
        sig = d[d.q < 0.05]
        m3 = sig[sig.n >= 3]
        summ[layer] = {
            "significant": int(len(sig)), "with_3plus_instruments": int(len(m3)),
            "egger_intercept_ok": int((m3.egger_intercept_p >= 0.05).sum()) if len(m3) else 0,
            "weighted_median_agrees": int(((np.sign(m3.wmedian) == np.sign(m3.beta)) & (m3.wmedian_p < 0.05)).sum()) if len(m3) else 0,
            "heterogeneity_p05": int((sig.het_p < 0.05).sum()) if "het_p" in sig else 0,
            "steiger_correct": int(sig.steiger_correct_direction.sum()),
        }
    (PROCESSED / "mr_sensitivity_summary.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    return summ


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
    d = pd.read_parquet(PROCESSED / "mr_sensitivity.parquet")
    print(d[d.gene.isin(["IL12B", "IL10", "CD274", "TNFRSF6B", "IRF5", "GSDMB", "NXPE1"])].round(3).to_string())


def load(store=None) -> int:
    """One mr_sensitivity record per significant estimate (gene x layer)."""
    from datetime import date
    from .graph import GraphStore, evidence_id
    store = store or GraphStore()
    d = pd.read_parquet(PROCESSED / "mr_sensitivity.parquet")
    d = d[d.q < 0.05]
    n = 0

    def f(v, nd=3):
        return None if v is None or pd.isna(v) else round(float(v), nd)
    with store.batch():
        for x in d.itertuples():
            node = store.resolve(x.gene)
            if not node or node["type"] != "Gene":
                continue
            rec = {"gene": x.gene, "layer": x.layer, "study": x.studyId, "n_instruments": int(x.n),
                   "heterogeneity_p": f(getattr(x, "het_p", None), 4),
                   "egger_slope": f(getattr(x, "egger_slope", None)), "egger_intercept_p": f(getattr(x, "egger_intercept_p", None), 4),
                   "weighted_median": f(getattr(x, "wmedian", None)), "weighted_median_p": f(getattr(x, "wmedian_p", None), 4),
                   "steiger_correct_direction": bool(x.steiger_correct_direction),
                   "note": "Egger and weighted median need >= 3 instruments; Steiger r2 approximated from z and N"}
            ev = store.add_evidence(evidence_id("mr_sensitivity", f"{x.gene}|{x.layer}"), "mr_sensitivity",
                                    "DATASET:MR_UKBPPP_UC" if x.layer == "plasma_protein" else "DATASET:MR_GUT_EQTL_UC",
                                    date.today().isoformat(), **rec)
            store.add_claim(node["id"], "causal_estimate_for",
                            "DATASET:MR_UKBPPP_UC" if x.layer == "plasma_protein" else "DATASET:MR_GUT_EQTL_UC", ev)
            n += 1
    return n
