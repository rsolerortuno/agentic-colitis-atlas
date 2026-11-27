"""Direction check of semantic.py on each kind of signed evidence (extended 2026-11-23)."""
from ucatlas.semantic import _direction

PATHWAY = {"type": "pathway_activity", "name": "TNFa", "cohorts_up": 8, "cohorts_down": 0}
RESPONSE = {"type": "treatment_response", "gene": "OSM",
            "arms": {"ustekinumab": {"diff_nr_minus_r": 0.65, "padj": 2.6e-05},
                     "golimumab": {"diff_nr_minus_r": -0.6, "padj": 0.99}}}
MR = {"type": "mendelian_randomization", "gene": "CD274", "significant": True, "beta": -0.52}
META = {"type": "expression_meta", "gene": "CXCR2", "significant": True, "log2fc": 2.8, "layer": "rna"}
LIT = {"type": "literature_relation", "subject_text": "AQP9", "object_text": "UC",
       "sentence": "AQP9 was predominantly expressed in neutrophils and upregulated in the colon tissues of mice with UC."}


def test_pathway_activity_direction():
    assert _direction("The TNFa pathway is decreased versus control in all eight cohorts.", None, [PATHWAY])
    assert not _direction("The TNFa pathway is increased versus control in all eight cohorts.", None, [PATHWAY])


def test_pathway_names_are_not_read_from_other_clauses():
    # "pathways rise together" is about pathways in general, not the one pathway named elsewhere in the sentence
    ev = [{"type": "pathway_activity", "name": "p53", "cohorts_up": 0, "cohorts_down": 3}]
    assert not _direction("So many pathways rise together, but p53 is not increased in any cohort.", None, ev)


def test_treatment_response_direction():
    assert _direction("OSM is lower at baseline in ustekinumab non-responders.", "OSM", [RESPONSE])
    assert _direction("OSM is higher in ustekinumab responders.", "OSM", [RESPONSE])
    assert not _direction("OSM is higher at baseline in ustekinumab non-responders.", "OSM", [RESPONSE])
    # a non-significant arm gives no direction
    assert not _direction("OSM is higher in golimumab non-responders.", "OSM", [RESPONSE])


def test_mendelian_randomisation_direction():
    assert _direction("Genetically lower plasma protein lowers UC risk.", "CD274", [MR])
    assert _direction("Genetically higher plasma PD-L1 raises UC risk.", "CD274", [MR])
    assert not _direction("Genetically lower plasma PD-L1 raises UC risk.", "CD274", [MR])
    assert not _direction("Genetically higher plasma PD-L1 is protective.", "CD274", [MR])
    # a risk sentence is not read as a disease-versus-control expression claim
    assert not _direction("Lower plasma protein raises UC risk.", "CD274", [MR, {**META, "gene": "CD274"}])


def test_ranking_word_exempts_only_its_own_clause():
    text = "CXCR2 tops the convergence ranking: pooled UC vs control log2FC is down in all eight cohorts."
    assert _direction(text, "CXCR2", [META])
    assert not _direction("CXCR2 is placed lower in the ranking.", "CXCR2", [META])


def test_score_elsewhere_in_the_sentence_does_not_hide_a_direction():
    text = "CXCR1 is a CXCR2 partner with expression down across cohorts and a moderate genetic score."
    assert _direction(text, "CXCR1", [{**META, "gene": "CXCR1"}])


def test_literature_sentence_direction():
    assert _direction("AQP9 is downregulated in UC mouse colon.", "AQP9", [LIT])
    assert not _direction("AQP9 is upregulated in UC mouse colon.", "AQP9", [LIT])


def test_literature_lists_of_genes_and_interventions():
    lst = {"type": "literature_relation", "subject_text": "SLC26A3", "object_text": "UC",
           "sentence": "The expression levels of ACSF2, NPY, SLC26A3 and BRINP3 were significantly lower in UC "
                       "patients, while the expression levels of CCL20 and LCN2 were higher."}
    assert _direction("SLC26A3 is reported higher in UC than in controls.", "SLC26A3", [lst])
    assert not _direction("SLC26A3 is reported lower in UC than in controls.", "SLC26A3", [lst])
    # an experiment that knocks a gene down says nothing about its expression in disease
    ko = {"type": "literature_relation", "subject_text": "CXCL8", "object_text": "UC",
          "sentence": "Knockdown of CXCL8 reduced the expression of inflammatory genes in colitis."}
    assert not _direction("CXCL8 is upregulated in UC.", "CXCL8", [ko])
    # disease abbreviations are not genes
    uc = {"type": "literature_relation", "subject_text": "UC", "object_text": "IBD",
          "sentence": "UC expression of the marker was increased."}
    assert not _direction("Expression is decreased in UC.", None, [uc])


def test_negation_after_the_direction_word_does_not_negate_it():
    lit = {"type": "literature_relation", "subject_text": "DUOX2", "object_text": "UC",
           "sentence": "Proteins most highly elevated in inflamed UC biopsies were HLA-DRA and DUOX2."}
    assert _direction("DUOX2 appears among the most reduced proteins in inflamed UC biopsies without an "
                      "established inhibitor.", None, [lit])
    assert not _direction("DUOX2 is not reduced in inflamed UC biopsies.", None, [lit])


def test_clause_about_several_genes_and_carried_subject():
    mr2 = [MR, {**MR, "gene": "CXCL5", "beta": -0.45}]
    text = "For CD274 and CXCL5 genetics runs opposite to an inhibitor strategy, since genetically lower protein is protective."
    assert _direction(text, None, mr2)
    assert not _direction(text.replace("lower protein", "higher protein"), None, mr2)
    # genes that disagree: nothing is attributed
    mixed = [MR, {**MR, "gene": "CXCL5", "beta": 0.45}]
    assert not _direction("CD274 and CXCL5: genetically lower protein is protective.", None, mixed)


def test_drug_names_need_a_cited_record(monkeypatch):
    import re
    from ucatlas import semantic
    monkeypatch.setattr(semantic, "_DRUGS", re.compile(r"(?<![\w-])(reparixin|navarixin)(?![\w-])", re.I))
    drug = {"type": "open_targets_drug", "gene": "CXCR2", "drug": "NAVARIXIN", "max_clinical_stage": "PHASE_2"}
    paper = {"type": "drug_literature", "gene": "CXCR2", "drug": "navarixin", "pmid": 1}
    assert semantic.uncited_drugs("Navarixin reached phase 2.", [drug, paper]) == []
    assert semantic.uncited_drugs("Reparixin and navarixin reached late-stage trials.", [drug, paper]) == ["reparixin"]
    # since 2026-11-27 a database listing alone is not enough: the drug needs a paper or a treatment arm
    assert semantic.uncited_drugs("Navarixin reached phase 2.", [drug]) == ["navarixin"]
    arm = {"type": "treatment_response", "gene": "IL12B", "drug": "navarixin"}
    assert semantic.uncited_drugs("Navarixin reached phase 2.", [arm]) == []
    # names inside annotations are not text
    assert semantic.uncited_drugs("A drug [drug=reparixin] exists.", []) == []


def test_a_drug_effect_is_not_the_gene_s_expression():
    """"The antagonist lowered blood neutrophil counts" says what a drug does, not that CXCR2 falls in UC."""
    text = "It is up in all eight cohorts, and the antagonist MK-7123 lowered blood neutrophil counts in a COPD trial."
    assert not _direction(text, "CXCR2", [META | {"gene": "CXCR2"}])
    assert _direction("CXCR2 is lowered in inflamed mucosa.", "CXCR2", [META | {"gene": "CXCR2"}])
