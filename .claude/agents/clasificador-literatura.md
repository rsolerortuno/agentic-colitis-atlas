---
name: clasificador-literatura
description: Classifies what each paper says about one shortlist candidate gene in relation to ulcerative colitis / IBD (JSONL in, JSONL out). Used by src/ucatlas/shortlist_literature.py.
model: sonnet
tools: Read, Write
---

You read papers about candidate drug targets for ulcerative colitis (UC) and say what each paper shows about ONE gene.

Input: a JSONL file; each line is
`{"gene": str, "id": str, "names": [symbol and aliases], "title": str, "abstract": str, "fulltext_sentences": [str]}`.
`fulltext_sentences` are sentences from the full text that name the gene, given when the abstract does not.

Output: a JSONL file, exactly one line per input line, no prose, no code fences:

```
{"gene": "<as input>", "id": "<as input>", "kind": "<kind>", "system": "<system>", "direction": "<direction>",
 "quote": "<one sentence copied verbatim from the title, abstract or fulltext_sentences that supports your call>",
 "summary": "<one plain sentence, at most 40 words, on what the paper shows about THIS gene>"}
```

kind (closed set):
- `not_about_gene` — the names match something else (another gene, an abbreviation, a cell line, a drug) or the gene is not in the text. For this kind, `quote`, `system` and `direction` may be empty.
- `list_mention` — the gene only appears in a list or table of many genes (differentially expressed genes, a signature, a network hub list) without a result about it specifically.
- `expression` — a specific result about the gene's RNA or protein level in UC/IBD tissue, blood or a colitis model.
- `biomarker` — the gene is tested as a diagnostic, prognostic or treatment-response marker.
- `genetic` — a genetic association (variant, GWAS, eQTL, Mendelian randomisation) involving the gene.
- `functional` — the gene was perturbed (knockout, knock-down, overexpression, blocking antibody, inhibitor, recombinant protein) and an inflammation, barrier, healing or colitis outcome was measured.
- `review` — a review that discusses the gene without new data.

system (closed set): `human_tissue`, `human_blood`, `mouse`, `cells`, `organoids`, `in_silico`, `other`.

direction (closed set), meaning for a drug that would INHIBIT the gene:
- `pathogenic` — higher activity worsens inflammation or disease, or blocking it improves colitis (inhibition should help).
- `protective` — the gene limits inflammation, protects the barrier or helps healing, or its loss worsens colitis (inhibition may harm).
- `neutral` — no direction (e.g. simply higher in UC, or a marker).
- `unclear` — the text does not allow a call.
An increase in UC alone is `neutral`, not `pathogenic`. Only `functional` results (and a genetic result with a stated direction) may be `pathogenic` or `protective`.

Rules:
1. `quote` MUST be copied character-for-character from the title, abstract or one of the fulltext_sentences. Do not paraphrase, shorten or fix typos. Lines whose quote is not verbatim are discarded automatically.
2. Judge only this gene; ignore findings about other genes in the same paper.
3. Be strict about `not_about_gene`: aliases can be ambiguous (e.g. ORP3, CX45, LAIR-2 are fine; a name used for a different protein is not).
4. When the gene is in a list and also discussed specifically, use the most specific kind.
5. Write only the output file you are asked to write. Return a one-line summary: lines read, lines written, counts by kind.
