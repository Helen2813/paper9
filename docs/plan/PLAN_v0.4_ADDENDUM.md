# Paper 9 — PLAN v0.4 addendum (2026-11-08): biological and general strengthening

Context: v030-v037 established the methodological core (the prioritized cell type depends on the
marker-gene selection rule; |log2FC| sets are mostly depleted genes; detection without resolution
among stromal states; Ochoa / tubal-epithelium claims not supported). Two gaps remain:
(1) little biology beyond "stromal, not resolvable"; (2) generality beyond one disease / atlas.

## Biology arm (locus level) — v039
Aggregate enrichment averages over loci; if risk loci act in different compartments, aggregate
cell-type rankings are unstable by construction. v039 maps each genome-wide significant locus
(KOLLER_EUR primary; GCST90205183, leiomyoma, hypertension, T2D for comparison) to an endometrial
compartment through its candidate genes (HECA), tests compartment enrichment against all expressed
genes, identifies the loci that drive each aggregate ranking (leave-one-locus-out), and flags
loci shared with leiomyoma / T2D / hypertension.
Possible next steps depending on v039: locus-level tubal (Weigert) and tissue (GTEx) context;
literature annotation of driver loci.

## Generality arm — v040 (planned)
(a) Positive-control traits with known cell types in the same atlas (hypertension -> vascular /
perivascular; leiomyoma -> smooth muscle): rank of the expected group under |log2FC|, up-only and
down-only sets. (b) Second atlas (Weigert fallopian tube): Koller-style marker sets per tube cell
type and the same rule comparison for the same traits.

No overlap with Paper 7 (progesterone-response timing) or Paper 8 (RNA markers, assayability).

## External review (2026-11-08) and response

Points accepted:
- Related work: Li et al. 2026 (medRxiv 10.1101/2025.05.24.25328275) benchmarked 20 GWAS-single-cell
  configurations with simulated GWAS (known beta-cell source), FWER under the null, signal
  redirection to related cell types when the causal type is missing, and the CATCH combination
  (conLDSC-Cepo + scDRS-mBAT-combo). They use top-10% specificity gene sets from eight metrics and
  do not examine directional (up / down / absolute) marker selection. Paper 9's contribution is
  framed as: how directional marker definition changes the biological reading of heritability
  enrichment, in a published case, and whether related endometrial states can be told apart.
  Working title: "Directional marker selection alters genetic cell-state prioritization: a
  reproducibility and calibration study in endometriosis".
- Koller cohort not yet exact (12-14 vs 22 endometriosis samples): K1b in v040 reconstructs it from
  the E-MTAB-14039 sample table before any claim about the published analysis.
- The planted-profile simulations may favour up-only sets by design: v041 adds truths independent of
  the marker rule (negative / depletion profiles, expression-level profiles).
- Calibration precision: v041 uses several hundred null rotations and reports FWER with 95%
  Clopper-Pearson intervals at 5% (not 10%); FDR per set and FWER per family kept separate.
- Gene-set size: v041 adds size-matched ABS / UP / DOWN sets (same N per type) and continuous
  specificity (no threshold).
- Ochoa 2025: wording stays "not reproduced under our reconstruction" (labels, cells and gene
  scoring differ from the original code); no refutation claim.
Planned order: v040 (locus atlas fix + K1b) -> v041 matched directional benchmark -> v042 external
validation (positive-control traits; second atlas) -> manuscript.
