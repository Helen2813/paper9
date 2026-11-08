# Paper 9 — PLAN v0.3 (frozen 2026-11-06)

Working title: Resolution limits of GWAS-based cell-state prioritization: a cross-atlas and
cross-trait audit in endometriosis.

Status of the v0.2 plan: all original hypotheses (G1-G6, X1-X4, H-sub) and their pre-registered
results stay in the registry and RESULTS_LOG unchanged and are reported as the original plan.
Everything below is post-review (labelled so in the paper).

## Core question (Option A, main)

Under which conditions can GWAS summary statistics reliably prioritize one cell state over
closely related ones, and how can unreliable prioritizations be detected?

Components (each fixed in CHANGELOG before it is run):
A1 Resolution: can decidual stroma be distinguished from perivascular / arterial / smooth-muscle
   programs? Covariate correlation and VIF; pairwise conditional tests; delete-one-block
   jackknife; simulations with a known source of signal (misassignment rates).
A2 Calibration: false-positive rate of parametric gene-property tests for donor-estimated
   contrasts under a true null (random splits within one tissue context); permutation
   calibration with enough permutations; jackknife-based cross-trait difference tests.
A3 Specificity: a small, pre-defined panel of control GWAS (hypertension, T2D, leiomyoma; more
   only if a systematic problem appears) run through every analysis.
A4 Published claims reproduced as closely as their methods allow: HECA 2024 (fGWAS-like),
   Koller 2026 (HECA marker genes |log2FC| > 1.2, FDR < 0.01, cases / controls separately,
   endometrial / immune separately, S-LDSC), Ochoa 2025 (scDRS, top 1,000 genes, lesion cells,
   their sample re-classification), then perturbed one factor at a time.
A5 Detection limits for every null.

## Exploratory branch (Option B, does not decide the paper)

B1 Are endometriosis risk associations preferentially enriched in fallopian-tube secretory or
   ciliated epithelial programs relative to matched endometrial epithelial programs?
   Direct difference / conditional tests, never rank comparisons; no claim about lesion origin.
   Reference: an independent fallopian-tube single-cell atlas (to be identified and fixed in
   CHANGELOG before any genetics is run on it).

## Dropped or reduced

G2 (TWAS) only if it answers a specific localization question. G6 stays exploratory context.

## Boundaries with Paper 7 / Paper 8

No PGR-response tests (Paper 7). No new central claims on differential expression, methylation,
lncRNA or the reserved RNA markers (Paper 8). Shared datasets serve only as references for
interpreting genetic signal.
