"""Plain-language text blocks of the report: tab intros and the one-line tab descriptions of the Sections menu.

Written for someone who knows a little biology but not this pipeline. Nothing here is computed.
"""

PLAIN: dict[str, str] = {
    "candidates": "This tab ranks every gene by how many separate lines of evidence point to it. A gene scores high when different kinds of data agree: it changes in UC biopsies, the change repeats across eight patient cohorts, it differs between inflamed and non-inflamed mucosa of the same patient, human genetics points to it, papers support it and a drug exists or looks possible. The table lists the top 25 and the plots show where each score comes from.",
    "expression": "This tab is the starting point of the analysis. It compares colon biopsies of people with ulcerative colitis with colon of people without inflammatory bowel disease, in eight public cohorts. A gene much higher in inflamed tissue may be a target or a marker. A gene much lower is usually a normal colon gene the inflamed lining has lost.",
    "cells": "A biopsy is a mix of lining cells, connective tissue, vessels and immune cells, and inflamed colon fills with neutrophils from the blood. A gene that looks high in the mix may come from any of them. This tab uses single-cell data to say which cell makes each gene, a composition adjustment to say how much of the change is cells arriving, and tissue sections to say where the gene sits and what is next to it.",
    "shortlist": "If you want to pick a new target to test in the laboratory, start here. Each first-in-class gene is checked against the questions you would ask before an experiment: is the increase reproducible and seen as protein, is it a change of the tissue, do we know which cell makes it, does genetics say it is a cause, can a drug reach it, is it safe, what do the papers on it say, and has anyone tried it in mice. Each gene also gets a suggested first experiment.",
    "known": "Some genes are clearly part of ulcerative colitis, with dozens or hundreds of papers, and still no drug reaches them. The first-in-class list leaves them out on purpose because it looks for new genes, and the repurposing list leaves them out because they have no drug at all. For someone who wants to make a drug, these may be the clearest openings: the biology is established and nobody owns the target yet. Each gene gets the same traffic lights as the other lists.",
    "repurposing": "The Validation shortlist looks for targets that no drug reaches yet. This page looks at the opposite case: genes that change in the inflamed colon and already have a drug that blocks them, approved or in trials for another disease, but never tried in ulcerative colitis. Testing an existing drug can be faster, but it needs the same proof, and a drug made for cancer may be too harsh for a lifelong disease. Each gene gets the same traffic lights as the shortlist and a suggested first test with that drug.",
    "viewer": "This is the one place where you can look at real tissue. Pick a gene and see where it is switched on across a colon section, then click a spot to see which cells are around it.",
    "pathways": "Genes act in pathways and cells talk to each other. This tab shows which signalling pathways and transcription factors are switched on in inflamed mucosa in all eight cohorts, and which cell sends which signal to which cell. Co-expression is not proof that a signal is used.",
    "genetics": "People inherit small DNA differences that raise or lower how much of a protein or an RNA they make, all their life. If people with the higher-level version get ulcerative colitis more often, the gene is more likely a cause than a consequence of the disease (Mendelian randomisation). This tab does that with protein levels in blood and with gene expression measured in the colon wall and in single gut cells.",
    "structure": "Explore the predicted 3D shapes of the candidate proteins (AlphaFold) and the cavities a small molecule could use. Grey is the protein and red marks the selected pocket. Pocket scores describe shape; they do not predict that a drug will bind.",
    "discovery": "Most top-ranked genes are already well known. This tab looks for the opposite: genes that rise in UC mucosa as a change of the tissue rather than cells flooding in, with little UC literature and no clinical drug, that look druggable and are low in vital organs. Genes that pass every filter are called first-in-class opportunities. They are hypotheses to test, not proven targets.",
    "response": "Two questions that matter to patients and doctors. Does a gene help tell who will respond to which drug? And does the answer change with the patient: children or adults, rectum or the rest of the colon, mild or severe disease? Biopsies taken before treatment in four trials and cohorts answer the first; eight cohorts and biopsy sites answer the second.",
    "biologics": "Most UC drugs are antibodies or other large molecules that bind a protein outside the cell. They need a protein on the cell surface or secreted, higher in inflamed mucosa, present in most patients and absent from vital organs. The first half of this tab scores genes on that. The second half asks which genes belong to the immune cells in the mucosa, since those are the classic targets in UC.",
    "validation": "A ranking is only worth trusting if it survives hard questions. This tab checks the candidates against other cohorts and other kinds of data, shows how much the order changes when the weights change, and tests whether the scores recover real clinical targets they were never shown.",
    "literature": "This tab covers what is already known. The literature part shows what PubMed abstracts say about the top genes. The drug part shows which drugs already exist, how far they got in clinical trials and how druggable each protein looks.",
    "network": "The network draws the top genes, the drugs that act on them, the genes known to interact and their link to ulcerative colitis. It helps you see clusters and neighbours. It is not proof of mechanism.",
    "agents": "The atlas is written by a team of AI agents rather than by one. This tab explains how that was tested: the same analysis run by one agent and by a team, how the reports were marked, what each set-up cost and what the checks found. It is about how the atlas is made, not about any gene.",
    "methods": "This tab holds the receipts: every dataset with its checksum, the models used, the limits of the analysis, and a searchable list of the evidence records the report cites.",
}

TOUR: list[tuple[str, str, str]] = [
    ("candidates", "Candidates", "Which genes have the most agreeing evidence? The ranking, how each score is built and how stable the order is."),
    ("known", "Well known in UC, no drug yet", "Genes with established UC biology that no drug reaches yet: the same lights as the shortlist."),
    ("shortlist", "First-in-class shortlist", "New targets with few papers and no drug: a traffic light per criterion and a proposed first experiment."),
    ("repurposing", "Repurposing", "Genes with an inhibitor already in the clinic for another disease, never tried in UC: the same lights."),
    ("expression", "Expression", "What changes in inflamed colon, how it repeats across eight cohorts and how it follows disease activity."),
    ("cells", "Cells & space", "Which cell makes each gene, how much of the change is neutrophils arriving, and what sits next to it in tissue."),
    ("viewer", "Tissue viewer", "Real colon sections with their H&E picture. Pick a gene, click a spot, see its neighbours."),
    ("pathways", "Pathways", "Immune pathways, transcription factors and cell-cell circuits in inflamed mucosa."),
    ("genetics", "Genetics", "Does genetically higher protein in blood, or expression in the colon, change UC risk?"),
    ("structure", "Structure", "Do the proteins have pockets a small molecule could use? A heuristic, shown with its check."),
    ("discovery", "Discovery", "Which genes are new, undrugged and still safe? The first-in-class screen with all its filters."),
    ("response", "Treatment & subgroups", "Do genes track response to four biologics, and does the answer change with age, site or activity?"),
    ("biologics", "Biologics & immune", "Which proteins suit an antibody, and which genes belong to the immune cells of the mucosa?"),
    ("validation", "Validation & robustness", "Do the calls hold in other cohorts and under different weights? Recovery of known UC targets."),
    ("literature", "Literature & drugs", "What do papers say, and which drugs already exist?"),
    ("network", "Network", "Genes, drugs and interactions drawn as a graph."),
    ("agents", "Agent benchmark", "Is a team of AI agents better than one? The experiment, the costs and the checks."),
    ("methods", "Methods & evidence", "Datasets, limits and every cited evidence record."),
]

NAV: list[tuple[str, list[tuple[str, list[str]]]]] = [
    ("Start", [("Start here", ["overview"])]),
    ("Findings", [("Ranking", ["candidates"]), ("New targets", ["discovery", "known", "shortlist", "repurposing"])]),
    ("Evidence", [("Molecular data", ["expression", "pathways"]), ("Cells and tissue", ["cells", "viewer"]),
                  ("Causal genetics", ["genetics"]), ("Protein structure", ["structure"])]),
    ("Clinical use", [("Patients and trials", ["response"]), ("Drug types", ["biologics"]),
                      ("What is already known", ["literature", "network"])]),
    ("Trust and methods", [("Checks", ["validation"]), ("How the atlas was made", ["agents", "methods"])]),
]
NAV_TEXT: dict[str, tuple[str, str]] = {"overview": ("Overview", "What the atlas is, how it was built and the answer in brief.")}
NAV_TEXT.update({k: (title, text) for k, title, text in TOUR})
