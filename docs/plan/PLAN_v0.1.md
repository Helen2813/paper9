# Paper 9 — Research Plan v0.1 (pre-preflight draft)

Status: DRAFT. Becomes FROZEN as v0.2 after preflight results are reviewed.
Any change after v0.2 is logged in `docs/plan/CHANGELOG.md` with reason and date,
before the corresponding analysis is run.

## Central question

Where and when does inherited endometriosis risk act — in the source tissue
(eutopic endometrium), in the lesion niche, or outside the uterus — and does this
differ between disease subtypes?

Rationale: germline genotype is fixed at conception, so genetic anchoring is not
confounded by menstrual phase, control definition, laboratory or platform — the
problems that limit cross-sectional eutopic transcriptomics (Paper 7, Paper 8).

## Scope boundaries (non-overlap)

- Paper 7: identifiability of delay vs attenuation in decidualization time courses.
  Paper 9 does not re-analyse time-course geometry.
- Paper 8: specificity and transportability of published (menstrual-effluent and tissue)
  markers, adenomyosis/other-pathology comparators, epithelial pseudobulk of GSE266265,
  methylation (X25–X27), new marker discovery (X16–X21).
  Paper 9 does not test markers, does not build classifiers, and does not estimate
  disease effects inside GSE266265.
- Shared datasets are referenced from the other projects' local folders, never copied
  into this repository.

## Hypotheses

| ID | Bold hypothesis | Primary data | Primary test | Falsifier / alternative result |
|----|-----------------|--------------|--------------|--------------------------------|
| G1 | Endometriosis SNP-heritability is more enriched in lesion-context cell programs (ectopic stroma, mesothelium, peritoneal immune) than in the same cell types in eutopic endometrium | GWAS GCST90205183; FinnGen replication; GSE213216 (Tan 2022) | MAGMA gene-property + S-LDSC on context-specific gene sets (cell type within tissue context vs same cell type elsewhere) | Eutopic contexts equally or more enriched → genetic support for a source-defect model |
| G2 | Genetically predicted expression (TWAS z) is sign-concordant with lesion-vs-eutopic differences but not with eutopic case-vs-control differences | S-PrediXcan GTEx v8 (uterus, ovary, vagina, others); GSE141549; GSE51981 (replication, referenced from Paper 8) | Genome-wide sign/rank concordance with LD-block-preserving permutation | Concordance with eutopic differences ≥ lesion differences → eutopic changes partly causal |
| G3 | Ovarian, peritoneal and deep endometriosis have partly distinct genetic architecture (rg < 1) and distinct cell-context enrichment | FinnGen subtype endpoints | LDSC rg between subtypes; G1 pipeline per subtype | rg CI includes 1 and enrichment profiles match |
| G4 | Risk enrichment concentrates in menstrual/shedding/early-proliferative cell states, not secretory/decidual states (genetic test of retrograde-menstruation timing) | Phase-resolved states (HECA; GSE266265 export from Paper 8); decidualization-response genes (GSE246591, referenced from Paper 7) | MAGMA/S-LDSC on phase-specific state gene sets | Decidual states equally enriched → genetic support for impaired decidualization as a cause |
| G5 | Cell-context enrichment profile is specific to endometriosis | FinnGen uterine leiomyoma; N80.0 (endometriosis of uterus, adenomyosis code); non-reproductive control trait | Same pipeline as G1/G4 | Profiles shared with leiomyoma → specificity claims in G1–G4 withdrawn |

## Fixed analysis rules (to be confirmed at v0.2)

- Primary GWAS: GCST90205183 (European-dominant meta-analysis excluding 23andMe).
  FinnGen is used for replication (G1, G2, G4) and as the only source for G3 and G5.
  Sample overlap between FinnGen and earlier meta-analyses is checked before any combined use.
- Gene definition: protein-coding genes, GRCh37/GRCh38 matched to the GWAS build;
  MAGMA window 35 kb upstream / 10 kb downstream (sensitivity: 0/0 and 100/100).
- Cell-context gene sets: top 10% specificity genes per cell type per context,
  minimum 50 cells per type per context and minimum 3 donors; donor-level pseudobulk only.
- Multiple testing: Holm within each hypothesis' primary family; secondary analyses
  reported with FDR and labelled exploratory.
- Negative controls: random gene sets matched on gene length and expression level;
  non-reproductive trait GWAS.
- Pivot rule: if no endometrial or lesion context is enriched in G1 and G4, the paper is
  reframed as "endometriosis genetic risk acts outside the endometrium" and extended to
  ovarian, immune and peritoneal contexts. This pivot is pre-registered here.

## Compute

MacBook Air (Apple Silicon). All analyses run on summary statistics and gene sets;
single-cell atlases are reduced to per-donor pseudobulk / specificity tables once.
No full atlas is loaded into memory.

## Version log

- v0.1 (2026-11-05): hypotheses G1–G5 drafted after reading Paper 7 preprint and
  Paper 8 working log v2.32–v2.64; earlier H1/H3/H4 dropped as Paper 8 territory.
