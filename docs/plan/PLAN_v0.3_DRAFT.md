> SUPERSEDED on 2026-11-06 by PLAN_v0.3.md (frozen). Kept for the record.

# Paper 9 — PLAN v0.3 (DRAFT, not frozen)

Status: draft written 2026-11-06 after two external reviews. Frozen only after the v021 audit.
The v0.2 hypotheses (G1-G5, X1-X3, H-sub) and their pre-registered results stay as they are and are
reported as the original plan; everything below that is new is labelled "post-review".

## Question

How robust is cell-context prioritization of inherited endometriosis risk across GWAS,
single-cell references and analytical frameworks?

## Design (Option A, post-review): a multiverse / specification analysis

Claims under test (published): HECA 2024 decidual stroma + macrophages; Ochoa 2025 lesion M2
macrophages (+ endothelial, DC); Guare 2025 proliferative stroma; Koller 2026 (if a cell-level
claim exists).

Factors varied, each level fixed before it is run:
1. GWAS: GCST90205183 (primary, public); FinnGen R13 ENDO_ALL (nested in GCST, not independent);
   subtypes (OVARY, PERITONEUM, DEEP); a larger GWAS if public summary statistics exist.
2. Reference / labels: HECA (authors' labels); GSE266265 (Tan-nomenclature transfer);
   GSE213216 with the authors' (Fonseca) labels; GSE213216 with our label transfer.
3. Context: eutopic only; lesion only; lesion minus eutopic (G1).
4. Framework: MAGMA gene-property (specificity | average); LDSC-SEG (top 10%); scDRS (top 1,000
   genes, default controls); fGWAS-like mean expression without conditioning (G4 check 5).
5. Background / conditioning: average expression; + constraint (LOEUF) + breadth (tau).

Read-outs per claim: share of specifications in which the claimed cell context is significant
after the pre-specified correction and has the claimed sign; the same share for unrelated traits
(HYPERTENSION, T2D) and the positive controls (LEIOMYOMA -> uSMCs, HYPERTENSION -> vascular);
detection limits (MDE) for every null.

A claim counts as robust only if it is significant in most specifications AND not reproduced by
unrelated traits AND positive controls pass in the same specifications.

## Bounded extension (Option B, post-review, exploratory)

A pre-specified list of GWAS loci (Rahmioglu 2023 / Koller 2026 lead loci): which cell contexts
express the prioritized genes, using published fine-mapping / colocalization only (no new causal
claims).

## Kept as exploratory

G3 subtype architecture (part B), X2 CC-GWAS, X3 ovarian tilt, G6 tissue level (GTEx).

## Dropped / re-scoped

G2 as "TWAS sign concordance" is re-scoped: concordance is descriptive, not causal; it stays
only if it adds to Option A (e.g., direction of prioritized genes in lesions vs eutopic).

## Boundaries with Paper 7 / Paper 8

No PGR-signature tests (Paper 7). No claims about marker DE, methylation change or lncRNA
(Paper 8); GSE141549 and GSE266265 are used only as context for genetic questions, with separate
code and targets.
