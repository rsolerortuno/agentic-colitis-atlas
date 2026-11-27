"""Semantic verification of statements against the evidence they cite (deterministic, no LLM).

Existence and gene linkage (verify.py) prove that a statement points at real evidence about the right gene. They do
not prove that the statement says what the evidence says. This module adds two checks:

1. Numbers. Match a specific metric in one unambiguous cited record, with its sign,
   scale and printed precision. Structured numeric_claims support nested metrics
   with explicit record context. Ranking/global number pools never supply support.
2. Direction. If a clause uses only "up" words (or only "down" words) and every signed evidence it cites for its
   subject points the other way, the statement contradicts its evidence. Signed evidence: disease-versus-control
   changes of genes, pathways and transcription factors; treatment non-responders versus responders; Mendelian
   randomisation effects on UC risk; the verbatim literature sentence cited (extended 2026-11-23, when the verifier
   sabotage test caught 26 % of flipped directions because only bulk expression was compared). Statements without a
   gene (summary, sections) are checked clause by clause for clauses that name exactly one subject.
"""
from __future__ import annotations

import math
import re


NUM = re.compile(r"(?<![\w.])([-−])?(\d+(?:\.\d+)?(?:[eE][-−+]?\d+)?)(\s?%)?(?!\w)")
# Words after which a number is a label, not a measurement. Whole words only: single letters ("p", "q") would exempt
# p- and q-values and any word ending in that letter ("from 0.5"). Numbers glued to letters ("M13", "CD70", "3p") are
# never extracted by NUM in the first place.
EXEMPT_BEFORE = re.compile(r"(\b(phase|rank|ranked|top|grade|stage|tier|version|step|chromosome|chr|pmid:?|uniprot|chembl)"
                           r"|#|[A-Za-z]-)\s*$", re.I)
THRESHOLD_BEFORE = re.compile(r"(\bat|<|≤|<=|>|≥|>=|threshold|cut-?off|alpha|α|fdr|bh)\s*$", re.I)
THRESHOLDS = {0.001, 0.01, 0.05, 0.1, 0.25}   # significance levels written as such ("not significant at 0.05")
DRUG_AFTER = re.compile(r"^\s+[A-Z][A-Z0-9-]{3,}")   # "IODINE I 124 GIRENTUXIMAB": part of a drug name
UP = re.compile(r"\b(up(?!\s+to\b)|upregulated|up-regulated|increased|elevated|higher|overexpress\w*|induced|gained|raised|rises?|rising)\b", re.I)
DOWN = re.compile(r"\b(down|downregulated|down-regulated|decreased|reduced|lower|lowered|depleted|suppressed|underexpress\w*|diminished)\b", re.I)
# A negated sentence ("is not higher") states the opposite of its direction word; it is skipped, not read literally.
NEGATION = re.compile(r"\b(not|no|never|neither|nor|without|n't)\b|n't\b", re.I)
# ordering words: "placed lower", "listed higher" describe a position in a list, not expression
DIRECTION_EXEMPT = re.compile(r"\b(survival|prognos\w*|hazard|dependency|rank\w*|score\w*|placed|listed|ordered|"
                              r"positioned|priorit\w*|moved|demoted|promoted)\b", re.I)


def numbers_in(text: str, *, with_spans: bool = False) -> list[dict]:
    out = []
    for m in NUM.finditer(text.replace("−", "-")):
        raw = m.group(2)
        before = text[max(0, m.start() - 12):m.start()]
        val = float(raw)
        decimals = len(raw.split(".")[1].split("e")[0].split("E")[0]) if "." in raw else 0
        sci = "e" in raw.lower()
        is_int = "." not in raw and not sci
        if is_int and 1900 <= val <= 2100 and not metric_field(before, text[m.end():m.end() + 60]):
            continue
        if EXEMPT_BEFORE.search(before) or DRUG_AFTER.match(text[m.end():m.end() + 20]):
            continue
        if val in THRESHOLDS and THRESHOLD_BEFORE.search(before):
            continue
        out.append({"text": m.group(0).strip(), "value": -val if m.group(1) else val, "decimals": decimals,
                    "sci": sci, "percent": bool(m.group(3))})
        if with_spans:
            out[-1].update(start=m.start(), end=m.end())
    return out


def _matches(n: dict, candidates) -> bool:
    v = n["value"]
    if not math.isfinite(v):
        return False
    for c in candidates:
        if not math.isfinite(c):
            continue
        target = c * 100 if n["percent"] else c
        if n["sci"]:
            # Scientific notation: tolerance follows the precision actually printed.
            mantissa = n["text"].lower().split("e")[0].lstrip("+-")
            digits = len(mantissa.split(".")[1]) if "." in mantissa else 0
            exponent = int(n["text"].lower().split("e")[1].strip().rstrip("%"))
            tol = 0.5 * 10 ** (exponent - digits)
        else:
            tol = 0.5 * 10 ** (-n["decimals"])
        if abs(v - target) <= tol + abs(target) * 1e-12:
            return True

    return False


# Only these source fields store fractions that can be displayed as percentages.
# A log fold-change, p-value or score must never silently become a percentage.
FRACTION_FIELDS = {"share_of_tumor_transcripts", "freq", "loss_freq", "gain_freq",
                   "p_top10", "p_top25", "p_top100", "cpm_share", "source_share"}
PERCENT_FIELDS = {"pct_expressing", "source_pct",
                  "tumour_cells_expressing_pct_weakest_quarter_of_donors",
                  "tumour_cells_expressing_pct_median_donor"}  # already stored on the 0–100 scale
CONTEXT_FIELDS = ("gene", "type", "source", "layer", "cohort", "cell_type", "celltype", "cell_group", "tissue", "arm", "contrast", "profile", "track")


def numeric_claim_text(claim: dict, evidence: dict) -> str:
    """Canonical, auditable annotation for a number, including its source context.

    Example claim: {evidence_id, field: ['arms', 'atezo', 'pvalue'],
                    display: '0.032', unit: 'native'}.
    The exact returned annotation must occur once in the statement. No free-form
    label can relabel a protein result as RNA, or a p-value as an effect size.
    """
    path = claim.get("field")
    if not isinstance(path, list) or not path or any(not isinstance(k, (str, int)) or isinstance(k, bool) for k in path):
        raise ValueError("field must be a nonempty path of keys/list indices")
    value = evidence
    try:
        for key in path:
            if isinstance(value, list) and (not isinstance(key, int) or key < 0):
                raise ValueError("invalid list index")
            value = value[key]
    except (KeyError, IndexError, TypeError) as exc:
        raise ValueError("metric path not found") from exc
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("metric is not a finite number")
    display = claim.get("display", "")
    if not isinstance(display, str) or not NUM.fullmatch(display):
        raise ValueError("display must be one numeric literal")
    unit = claim.get("unit")
    if unit not in ("native", "percent", "native_percent") or "%" in display:
        raise ValueError("unit must be native, percent or native_percent; display must omit %")
    if unit == "percent" and path[-1] not in FRACTION_FIELDS:
        raise ValueError("metric is not registered as a fraction")
    if unit == "native_percent" and path[-1] not in PERCENT_FIELDS:
        raise ValueError("metric is not registered as a native percentage")
    # Prefix makes all magnitudes, including small integers and years, measurable.
    tokens = numbers_in("log2FC " + display + ("%" if unit == "percent" else ""))
    if len(tokens) != 1 or not _matches(tokens[0], [value]):
        raise ValueError("display does not match metric sign/scale/precision")
    context = "; ".join(f"{k}={evidence[k]}" for k in CONTEXT_FIELDS if evidence.get(k) is not None)
    return f"[{context}; {'.'.join(map(str, path))}={display}{'%' if unit in ('percent', 'native_percent') else ''}]"


def check(text: str, gene: str | None, evidence: list[dict], extra_values: set[float] | None = None,
          numeric_claims: list[dict] | None = None) -> dict:
    """Return unmatched numbers and direction contradictions for one statement."""
    # Never accept a number merely because it occurs elsewhere in the evidence,
    # a ranking row, a date, an identifier or a global summary.
    unmatched, binding_errors = [], []
    remaining = text
    indexed = {e.get("evidence_id"): e for e in evidence if e.get("evidence_id")}
    if numeric_claims is not None and not isinstance(numeric_claims, list):
        binding_errors.append("numeric_claims must be a list")
        numeric_claims = []
    for i, claim in enumerate(numeric_claims or []):
        try:
            if not isinstance(claim, dict):
                raise ValueError("claim must be an object")
            ev = indexed.get(claim.get("evidence_id"))
            if ev is None or (gene and (ev.get("gene") or "").upper() != gene.upper()):
                raise ValueError("claim must reference cited evidence for the same gene")
            annotation = numeric_claim_text(claim, ev)
            if remaining.count(annotation) != 1:
                raise ValueError("canonical metric/context annotation must occur exactly once")
            remaining = remaining.replace(annotation, " " * len(annotation), 1)
        except (ValueError, TypeError) as exc:
            binding_errors.append(f"numeric_claims[{i}]: {exc}")
    for n in numbers_in(remaining, with_spans=True):
        start, end = n["start"], n["end"]
        before, after = remaining[max(0, start - 80):start], remaining[end:end + 60]
        field = metric_field(before, after)
        candidates = []
        context_mismatch = False
        layers = re.findall(r"\b(RNA|protein|proteome|phospho|phosphoproteome)\b", before, re.I)
        named_layer = {"proteome":"protein", "phosphoproteome":"phospho"}.get(layers[-1].lower(), layers[-1].lower()) if layers else None
        for e in evidence:
            if gene and (e.get("gene") or "").upper() != gene.upper():
                continue
            if field and isinstance(e.get(field), (int, float)) and not isinstance(e.get(field), bool):
                candidates.append(e[field])
                if named_layer and e.get("layer") and e["layer"] != named_layer:
                    context_mismatch = True
        if context_mismatch or len(candidates) != 1 or (n["percent"] and field not in FRACTION_FIELDS) or not _matches(n, candidates):
            unmatched.append(n["text"])
    return {"unmatched_numbers": unmatched, "direction_contradiction": _direction(text, gene, evidence),
            "numeric_binding_errors": binding_errors, "uncited_drugs": uncited_drugs(text, evidence)}


# ------------------------------------------------------------------ drug names (2026-11-23)
# A reading check of run 2026-11-14 found about 30 drug names and stages written from memory or from another tool's
# output without citing the record that holds them. A drug named in a statement must appear in a cited record (Open
# Targets drugs, clinical indications, a treatment arm, a literature sentence).
_DRUGS: re.Pattern | None = None


def _drug_pattern() -> re.Pattern | None:
    global _DRUGS
    if _DRUGS is None:
        import json as _json
        names: set[str] = set()
        try:
            from .graph import GraphStore
            conn = GraphStore().conn
            for (p,) in conn.execute("select props from evidence where type in ('open_targets_drug', 'clinical_indications')"):
                d = _json.loads(p)
                vals = [d.get("drug")] + [x for k, v in d.items() if k.endswith("_drugs") and isinstance(v, list) for x in v]
                names |= {str(v).strip().lower() for v in vals if v}
        except Exception:   # no graph (CI without data): nothing to check against
            names = set()
        names = {n for n in names if re.fullmatch(r"[a-z][a-z0-9-]{4,}", n)}
        _DRUGS = re.compile(r"(?<![\w-])(" + "|".join(sorted(map(re.escape, names), key=len, reverse=True)) + r")(?![\w-])",
                            re.I) if names else re.compile(r"(?!x)x")
    return _DRUGS


def uncited_drugs(text: str, evidence: list[dict]) -> list[str]:
    """Drug names in the text (outside annotations) that no cited record contains."""
    import json as _json
    found = {m.group(1).lower() for m in _drug_pattern().finditer(re.sub(r"\[[^\]]*\]", " ", text))}
    if not found:
        return []
    blob = _json.dumps(evidence, default=str).lower()
    return sorted(n for n in found if n not in blob)


SIG = 0.05


def _sign(x) -> int | None:
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x == 0:
        return None
    return 1 if x > 0 else -1


def _add(d: dict, subject: str, key: str, sign: int | None) -> None:
    if sign is None or not subject:
        return
    g = d.setdefault(subject.upper(), {})
    g.setdefault("*", set()).add(sign)
    g.setdefault(key, set()).add(sign)


def _de_signs(evidence: list[dict]) -> dict[str, dict[str, set[int]]]:
    """subject -> layer -> signs of the significant disease-versus-control changes cited (layer "*" pools all).

    Subjects are genes (bulk differential expression, the cross-cohort meta-analysis, the covariate- and
    composition-adjusted changes, the discovery profile) and pathways or transcription factors (activity in the
    eight cohorts, signed only when every cohort that moved moved the same way).
    """
    signs: dict[str, dict[str, set[int]]] = {}
    for e in evidence:
        t, gene = e.get("type"), e.get("gene")
        if t in ("omics_de", "expression_meta") and e.get("significant") and gene:
            _add(signs, gene, str(e.get("layer") or "*").lower(), _sign(e.get("log2fc")))
        elif t in ("expression_adjusted", "composition_effect") and gene:
            for pre in ("adjusted", "unadjusted"):
                p = e.get(f"{pre}_padj")
                if isinstance(p, (int, float)) and p < SIG:
                    _add(signs, gene, str(e.get("layer") or "*").lower(), _sign(e.get(f"{pre}_log2fc")))
        elif t == "discovery_profile" and gene and (e.get("cohorts_up") or 0) >= 5:
            _add(signs, gene, "rna", _sign(e.get("meta_log2fc")))
        elif t in ("pathway_activity", "tf_activity") and e.get("name"):
            up, down = e.get("cohorts_up") or 0, e.get("cohorts_down") or 0
            if up != down and min(up, down) == 0:
                _add(signs, e["name"], "*", 1 if up else -1)
    return signs


def _response_signs(evidence: list[dict]) -> dict[str, dict[str, int]]:
    """subject -> treatment arm -> sign of non-responders minus responders at baseline (significant arms only)."""
    out: dict[str, dict[str, int]] = {}
    for e in evidence:
        if e.get("type") == "treatment_response" and e.get("gene"):
            subject, field = e["gene"], "diff_nr_minus_r"
        elif e.get("type") == "response_pathway" and e.get("name"):
            subject, field = e["name"], "t_nr_vs_r"
        else:
            continue
        for arm, a in (e.get("arms") or {}).items():
            s = _sign(a.get(field)) if isinstance(a, dict) else None
            if s is not None and isinstance(a.get("padj"), (int, float)) and a["padj"] < SIG:
                out.setdefault(subject.upper(), {})[arm] = s
    return out


def _mr_signs(evidence: list[dict]) -> dict[str, dict[str, set[int]]]:
    """gene -> exposure layer (plasma, gut, "*") -> sign of the significant MR effect (beta > 0: higher raises risk)."""
    out: dict[str, dict[str, set[int]]] = {}
    for e in evidence:
        if e.get("type") in ("mendelian_randomization", "mendelian_randomization_gut") and e.get("significant") \
                and e.get("gene"):
            _add(out, e["gene"], "plasma" if e["type"] == "mendelian_randomization" else "gut", _sign(e.get("beta")))
    return out


EXPRESSION_CUE = re.compile(r"\b(express\w*|levels?|mRNA|transcripts?|proteins?|abundan\w*|upregulat\w*|"
                            r"downregulat\w*|up-regulat\w*|down-regulat\w*|overexpress\w*|underexpress\w*)\b", re.I)
INTERVENTION = re.compile(r"\b(knock\w*|deficien\w*|deplet\w*|inhibit\w*|silenc\w*|block\w*|antagon\w*|treat\w*|"
                          r"administ\w*|agonist\w*|transfect\w*|null|-/-|overexpression\s+of|si\w*RNA|shRNA)\b", re.I)


def _lit_signs(evidence: list[dict]) -> dict[str, set[int]]:
    """gene -> directions that the cited verbatim literature sentences give for the gene's expression in disease.

    A source clause counts only when it names the gene, speaks of expression, is not negated and does not describe
    an experimental intervention ("knockdown of X reduced ..."): what the paper measured, not what it did.
    """
    out: dict[str, set[int]] = {}
    for e in evidence:
        if e.get("type") != "literature_relation" or not e.get("sentence"):
            continue
        for name in {str(e.get(k) or "") for k in ("subject_text", "object_text", "gene")}:
            if not re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", name) or name in NOT_GENES:
                continue                       # gene symbols only, not "IBD" or "UC"
            # source sentences list genes ("ACSF2, NPY, SLC26A3 and PKLPP2 were lower"): split at clause breaks only
            for clause in LIT_CLAUSE.split(e["sentence"]):
                if not re.search(rf"(?<![\w-]){re.escape(name)}(?![\w-])", clause) or not EXPRESSION_CUE.search(clause):
                    continue
                if NEGATION.search(clause) or INTERVENTION.search(clause) or DIRECTION_EXEMPT.search(clause):
                    continue
                said = _said(clause)
                if said is not None:
                    out.setdefault(name.upper(), set()).add(said)
    return out


# Clauses are checked one at a time, so "protein down ... but RNA up" or a survival clause next to an expression clause
# no longer hides the whole sentence (verifier sabotage test, 2026-11-03: 16 % of flipped directions were caught).
CLAUSE = re.compile(r"[;,:]|\bbut\b|\bwhile\b|\bwhereas\b|\band\b|\bwith\b", re.I)
LIT_CLAUSE = re.compile(r"[;:]|\bbut\b|\bwhile\b|\bwhereas\b", re.I)
NOT_GENES = {"UC", "IBD", "CD", "HC", "NC", "DSS", "TNBS", "RNA", "DNA", "AUC", "ROC"}
# Survival or dependency words make a whole sentence prognostic; ordering words ("ranks first", "placed above",
# "score") only exempt their own clause ("tops the ranking: log2FC is up" still states a direction).
SENTENCE_EXEMPT = re.compile(r"\b(survival|prognos\w*|hazard|dependency)\b", re.I)
LAYER_WORD = re.compile(r"\b(RNA|mRNA|transcripts?|protein|proteome|phospho\w*)\b", re.I)
PRINTED_LFC = re.compile(r"log2\s*FC\s*[=:]?\s*([-−]?\d+(?:\.\d+)?)", re.I)
# "control" side of a disease-vs-control comparison (healthy mucosa, non-inflamed mucosa, adjacent normal tissue)
CONTROL_WORDS = r"(normal|healthy|adjacent|controls?|non-?inflamed|uninflamed|quiescent)"
RELATIVE_TO_NORMAL = re.compile(r"\b(higher|lower|up|down|increased|decreased|elevated|reduced)\s+in\s+(the\s+)?"
                                + CONTROL_WORDS, re.I)
NORMAL_ONLY = re.compile(r"\b" + CONTROL_WORDS + r"\b", re.I)
# "decreased versus control": the control is the reference, so the clause is about disease
VERSUS_CONTROL = re.compile(r"\b(versus|vs\.?|compared\s+(with|to)|relative\s+to|than)\s+(the\s+)?(healthy\s+|normal\s+)?"
                            + CONTROL_WORDS + r"\b", re.I)
# "disease" side: tumour (oncology) or inflamed / lesional / diseased tissue (immune-mediated disease)
TUMOUR = re.compile(r"\b(tumou?rs?|(?<!non-)(?<!non)inflamed|inflammation|lesions?|lesional|diseased|active disease|"
                    r"colitis|UC)\b", re.I)
# treatment response: the side the direction word refers to ("lower in non-responders", "higher in responders")
RESPONDER = re.compile(r"\b(non-?responders?|responders?)\b", re.I)
# genetic risk: Mendelian randomisation sentences ("genetically lower protein raises UC risk", "... is protective")
RISK = re.compile(r"\b(risk|protective|protects?|susceptibility)\b", re.I)
OUT_UP = re.compile(r"\b(raises|raising|increases|increasing)\b|\b(higher|increased|greater|elevated)\s+(\w+\s+)?risk\b"
                    r"|\brisk-increasing\b|\bharmful\b", re.I)
OUT_DOWN = re.compile(r"\b(lowers|lowering|reduces|reducing|decreases|decreasing|protect\w*)\b"
                      r"|\b(lower|reduced|decreased|smaller)\s+(\w+\s+)?risk\b|\brisk-reducing\b", re.I)


def _said(clause: str) -> int | None:
    """+1 / -1 when the clause has only "up" or only "down" words, read against healthy tissue; None otherwise."""
    up, down = bool(UP.search(clause)), bool(DOWN.search(clause))
    if up == down:
        return None
    said = 1 if up else -1
    if RELATIVE_TO_NORMAL.search(clause):
        said = -said                      # "lower in healthy mucosa" = higher in disease
    elif NORMAL_ONLY.search(clause) and not TUMOUR.search(clause) and not VERSUS_CONTROL.search(clause):
        return None
    return said


def _layer_of(clause: str) -> str | None:
    m = LAYER_WORD.findall(clause)
    if not m:
        return None
    w = m[-1].lower()
    return "rna" if w in ("rna", "mrna", "transcript", "transcripts") else "protein" if w.startswith("prot") else "phospho"


# a clause about what a drug does ("the antagonist lowered blood neutrophil counts") is not a statement about the
# gene's expression in disease (2026-11-26)
DRUG_EFFECT = re.compile(r"\b(antagonists?|inhibitors?|blockers?|antibod(y|ies))\b", re.I)
DIRECTION_WORD = re.compile(UP.pattern + "|" + DOWN.pattern + "|" + OUT_UP.pattern + "|" + OUT_DOWN.pattern, re.I)


def _drug_effect(clause: str) -> bool:
    """The drug is the subject that acts ("the antagonist lowered ...") when it comes before the direction word;
    "most reduced ... without an inhibitor" stays a statement about the gene."""
    m, d = DRUG_EFFECT.search(clause), DIRECTION_WORD.search(clause)
    return bool(m and d and m.start() < d.start())


def _negated(clause: str) -> bool:
    """A negation before the first direction word ("is not higher"); one after it ("higher, without an inhibitor")
    does not negate the direction."""
    n, d = NEGATION.search(clause), DIRECTION_WORD.search(clause)
    return bool(n) and (d is None or n.start() < d.start())


# Each clause check returns (verdict, message): verdict "contra" (the clause contradicts the evidence), "ok" (it agrees)
# or None (nothing to compare).
def _response_clause(clause: str, arms: dict[str, int]) -> tuple[str | None, str]:
    m = re.search(r"\b(higher|lower|up|down|increased|decreased|elevated|reduced)\b.*?" + RESPONDER.pattern, clause, re.I)
    if not m or _negated(clause) or not arms:
        return None, ""
    said = 1 if UP.fullmatch(m.group(1)) else -1
    if not re.match(r"non", m.group(2), re.I):
        said = -said                      # "higher in responders" = lower in non-responders
    named = [a for a in arms if re.search(rf"\b{re.escape(a.split(' (')[0])}\b", clause, re.I)]
    seen = {arms[a] for a in (named or arms)}
    if len(seen) != 1:
        return None, ""
    if said != next(iter(seen)):
        side = "non-responders" if said > 0 else "responders"
        return "contra", f"text says higher in {side}, cited treatment-response evidence says the opposite"
    return "ok", ""


def _risk_clause(clause: str, mr: dict[str, set[int]]) -> tuple[str | None, str]:
    if _negated(clause) or not mr:
        return None, ""
    o_up, o_down = OUT_UP.search(clause), OUT_DOWN.search(clause)
    if bool(o_up) == bool(o_down):
        return None, ""
    rest = (OUT_UP if o_up else OUT_DOWN).sub(" ", clause)
    expo = _said(rest) if not RELATIVE_TO_NORMAL.search(rest) else None
    if expo is None:
        return None, ""
    implied = expo * (1 if o_up else -1)  # sign of the MR effect the sentence states
    layer = "plasma" if re.search(r"\b(plasma|circulating|serum|pQTL)\b", clause, re.I) else \
        "gut" if re.search(r"\b(gut|colon\w*|intestin\w*|mucos\w*|cell|eQTL|tissue)\b", clause, re.I) else "*"
    seen = mr.get(layer, mr.get("*", set()))
    if len(seen) != 1:
        return None, ""
    if implied != next(iter(seen)):
        return "contra", (f"text implies that higher protein or expression {'raises' if implied > 0 else 'lowers'} "
                          f"UC risk, cited Mendelian randomisation says it {'lowers' if implied > 0 else 'raises'} it")
    return "ok", ""


def _de_clause(clause: str, sentence: str, de: dict, lit: set[int]) -> tuple[str | None, str]:
    # A survival statement such as "higher expression predicts better survival"
    # is not a claim that expression increased in tumour versus normal tissue.
    # Keep this a narrow DE wording check, not a bag of direction words.
    if DIRECTION_EXEMPT.search(clause) or _negated(clause) or RISK.search(clause) or _drug_effect(clause):
        return None, ""
    # in a survival sentence only clauses that name the tumour, a layer or a log2FC are about
    # tumour-versus-normal expression ("higher expression is associated with worse survival" is not)
    if SENTENCE_EXEMPT.search(sentence) and not (TUMOUR.search(clause) or LAYER_WORD.search(clause)
                                                  or PRINTED_LFC.search(clause)):
        return None, ""
    said = _said(clause)
    if said is None:
        return None, ""
    word = "up" if said > 0 else "down"
    m = PRINTED_LFC.search(clause)
    if m and float(m.group(1).replace("−", "-")) != 0 and (float(m.group(1).replace("−", "-")) > 0) != (said > 0):
        return "contra", f"text says {word} next to a printed log2FC of the opposite sign ({m.group(1)})"
    layer = _layer_of(clause)
    seen_set = de.get(layer) if layer and layer in de else de.get("*", set())
    if len(seen_set) == 1:
        seen = next(iter(seen_set))
        return ("contra", f"text says {word}, cited omics evidence says {'up' if seen > 0 else 'down'}") \
            if said != seen else ("ok", "")
    if not seen_set and len(lit) == 1:
        return ("contra", f"text says {word}, the cited literature sentence says {'down' if said > 0 else 'up'}") \
            if said != next(iter(lit)) else ("ok", "")
    return None, ""


def _direction(text: str, gene: str | None, evidence: list[dict]) -> str | None:
    """Direction words that contradict the signed evidence cited.

    Four kinds of signed evidence: disease-versus-control changes of genes, pathways and transcription factors;
    baseline differences between treatment non-responders and responders; Mendelian randomisation effects on UC risk;
    and the direction of expression stated in the verbatim literature sentences cited. A clause about treatment
    response is read only against the response evidence, a clause about genetic risk only against Mendelian
    randomisation, any other clause against the disease-versus-control changes and the literature.

    With a statement gene every sentence is read as being about it. Without one (summary and section paragraphs) a
    clause is read as being about the subjects it names; a clause that names none is about the one gene the previous
    clause of its sentence named, or else the one gene its sentence names. A clause about several subjects is flagged
    only when it contradicts the evidence of every subject it names that has any ("for CD274 and CXCL5 higher protein
    is protective"), so a multi-gene comparison is never attributed to the wrong gene.
    """
    de, resp, mr, lit = _de_signs(evidence), _response_signs(evidence), _mr_signs(evidence), _lit_signs(evidence)
    # pathway and transcription-factor names must be named in the clause itself ("pathways rise together" is not p53)
    pathways = {str(e.get("name")).upper() for e in evidence if e.get("type") in ("pathway_activity", "tf_activity", "response_pathway")}
    subjects = set(de) | set(resp) | set(mr) | set(lit)
    if not subjects:
        return None

    def named_in(s: str) -> list[str]:
        return [g for g in subjects if re.search(rf"(?<![\w-]){re.escape(g)}(?![\w-])", s, re.I)]

    def check(clause: str, sentence: str, target: str) -> tuple[str | None, str]:
        if RESPONDER.search(clause):
            return _response_clause(clause, resp.get(target, {}))
        if RISK.search(clause) and target in mr:
            return _risk_clause(clause, mr[target])
        return _de_clause(clause, sentence, de.get(target, {}), lit.get(target, set()))

    for sentence in re.split(r"(?<=[.!?])\s+", text):
        in_sentence = [g for g in named_in(sentence) if g not in pathways]
        single = in_sentence[0] if len(in_sentence) == 1 else None
        last = None
        # canonical annotations use ";" between fields: keep each annotation inside its clause
        flat = re.sub(r"\[[^\]]*\]", lambda m: re.sub(r"[;,:]", " ", m.group(0)), sentence)
        for clause in CLAUSE.split(flat):
            if gene:
                targets = [gene.upper()]
            else:
                here = named_in(re.sub(r"\[[^\]]*\]", " ", clause))
                targets = here or ([last] if last else [single] if single else [])
                if len(here) == 1:
                    last = here[0] if here[0] not in pathways else None
                elif len(here) > 1:
                    last = None
            results = [(t, *check(clause, sentence, t)) for t in targets]
            contra = [(t, msg) for t, v, msg in results if v == "contra"]
            if contra and not any(v == "ok" for _, v, _ in results):
                t, msg = contra[0]
                where = "" if gene else f" ({', '.join(x for x, _ in contra)})"
                return msg + where
    return None


# Narrow grammar for legacy prose. Complex/nested comparisons must be emitted as
# structured claims rather than silently checked against a bag of numbers.
def metric_field(before: str, after: str) -> str | None:
    patterns = [
        (r"log2\s*FC\s*[=:]?\s*$", "log2fc"),
        (r"(?:padj|q(?:-value)?)\s*[=:]?\s*$", "padj"),
        (r"(?:pvalue|p-value|\bp)\s*[=:]?\s*$", "pvalue"),
        (r"(?:Cox\s+)?log\s*HR\s*[=:]?\s*$", "log_hr_per_sd"),
    ]
    for pattern, field in patterns:
        if re.search(pattern, before, re.I):
            return field
    if re.match(r"\s+pairs\b", after, re.I):
        return "n_pairs"
    if re.match(r"\s+of (?:its |the )?tumou?r transcripts", after, re.I):
        return "share_of_tumor_transcripts"
    return None
