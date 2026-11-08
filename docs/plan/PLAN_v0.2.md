# Paper 9 — Research Plan v0.2 (frozen after preflight v001, before any analysis)

Supersedes v0.1. Changes are listed in `CHANGELOG.md`. Any later change is logged there
with reason and date before the affected analysis is run.

## Central question

Where and when does inherited endometriosis risk act — in the source tissue
(eutopic endometrium), in the lesion niche, or outside the uterus — and does this
differ between disease subtypes?

## Prior result that constrains the plan

HECA (Marečková et al., Nat Genet 2024) ran fGWAS on an endometriosis GWAS meta-analysis and
reported enrichment in decidualized stromal cells (dStromal early/mid) and uterine
macrophages, using eutopic endometrium only (47 biopsies; 2 menstrual). GBMI (Guare et al.)
prioritised 35 eutopic cellular contexts. Therefore eutopic cell-type enrichment alone is
not new; Paper 9 is new where it (a) contrasts eutopic vs lesion vs peritoneal contexts,
(b) stress-tests the decidual-stroma result, (c) resolves subtypes, (d) tests sign
concordance of genetically predicted expression.

## Data (verified by preflight v001, 2026-11-05)

GWAS
- GCST90205183 (Rahmioglu 2023, excl. 23andMe): GRCh37 file + harmonised GRCh38 file;
  columns chr, pos, alleles, EAF, beta, SE, p (no per-SNP N; total sample size 474,160
  in metadata; case/control split taken from GWAS Catalog API in v002).
- FinnGen R13 (GRCh38, rsIDs; controls exclude other female genital disease for N14 codes):

| Endpoint | Cases | Controls | Use |
|---|---|---|---|
| N14_ENDOMETRIOSIS | 20,534 | 122,622 | all endometriosis |
| N14_ENDOMETRIOSIS_OVARY | 8,002 | 122,622 | G3 |
| N14_ENDOMETRIOSIS_PELVICPERITONEUM | 7,758 | 122,622 | G3 |
| N14_ENDOMETRIOSIS_DEEP | 3,890 | 271,625 | G3 |
| N14_ENDOMETRIOSIS_RECTPVAGSEPT_VAGINA | 3,263 | 122,622 | G3 (power-gated) |
| N14_ENDOMETRIOSIS_ASRM_STAGE3_4 | 10,284 | 271,625 | G3 severity |
| N14_ENDOMETRIOSIS_UTERUS (adenomyosis) | 5,864 | 122,622 | G5 |
| CD2_BENIGN_LEIOMYOMA_UTERI | 42,816 | 239,093 | G5 |
| I9_HYPTENS, T2D | — | — | G5 non-reproductive controls |

  Excluded: intestine (650), fallopian tube (284), skin scar (125), NOS (heterogeneous),
  C3_OVARY_ENDOMETROID (208; ovarian cancer, not an endometriosis subtype; too small).

Expression / atlases
- GSE213216 (Tan 2022): 51 samples, 21 patients; peritoneal endometriosis 25, eutopic 10,
  endometrioma 8, unaffected ovary 4, endometriosis-free peritoneum 4; RAW.tar 15.7 GB.
  Cell annotations: source to be confirmed in v002 (GEO RAW, HECA projection, or own annotation).
- GSE141549 (Gabriel 2020; local copy in Paper 8): GPL13376 n=336 (endometrium 105 = 62
  patient + 43 control; peritoneal lesion 72; DIE 33; endometrioma 28; sacrouterine 24;
  rectovaginal 22; peritoneum 52) and GPL10558 n=72. Cycle phase includes
  "medication" (140/336) — hormonal treatment is a recorded covariate.
- GSE266265 (60 donors; 35 ENDO / 25 non-ENDO): phase is not in GEO metadata but is in the
  RDS metadata exported by Paper 8 v2.32c (linked read-only).
- GSE51981 (148, GPL570; phases PE/ESE/MSE/LSE) — replication only, linked from Paper 8.
- GSE246591 decidualization time course — Paper 7 sample manifest linked; expression
  re-downloaded only if G4 needs it.

## Sample overlap (new, binding)

FinnGen contributed to the Rahmioglu 2023 meta-analysis (FinnGen Endometriosis Taskforce).
FinnGen R13 is therefore NOT an independent replication of GCST90205183. Rules:
- GCST90205183 is primary for G1, G2 and G4 (largest power).
- FinnGen R13 is primary for G3 and G5 (subtypes and comparator traits from one pipeline).
- Agreement between the two is reported as consistency, never as replication; the
  bivariate LDSC intercept quantifies overlap.

## Hypotheses (v0.2)

| ID | Hypothesis | Primary data | Primary test | Falsifier / alternative |
|----|-----------|--------------|--------------|-------------------------|
| G1 (main) | Heritability is more enriched in lesion-context cell programs (ectopic stroma, mesothelium, peritoneal immune) than in the same cell types in eutopic endometrium | GCST90205183; GSE213216 | S-LDSC (LDSC-SEG) and MAGMA gene-property on context-specific top-10% gene sets; contrast = lesion-context minus eutopic-context coefficient per cell type | Eutopic ≥ lesion → source-defect model supported |
| G2 | TWAS z (S-PrediXcan, GTEx v8 mashr: uterus, ovary, vagina + others) is sign-concordant with lesion-vs-own-eutopic differences but not with patient-vs-control eutopic differences | GCST90205183 (harmonised GRCh38); GSE141549 GPL13376 | Genome-wide Spearman / sign concordance; LD-block permutation; paired within patient (Sample_link); medication as covariate | Eutopic concordance ≥ lesion concordance |
| G3 | Ovarian, peritoneal and deep disease have partly distinct genetic architecture (rg < 1) and distinct cell-context profiles | FinnGen R13 subtypes | LDSC rg (Wald test of rg = 1); G1 pipeline per subtype; gate: h2 z ≥ 4 | rg CI includes 1 and same profiles |
| G4 (revised) | The published decidual-stroma enrichment (HECA) is not robust: it does not survive an independent phase-resolved atlas, a different method, removal of the top locus, or the lesion-stroma contrast | GCST90205183; GSE266265 non-ENDO donors (via Paper 8 export); HECA if downloadable | S-LDSC/MAGMA on phase × cell-type sets; leave-one-locus-out (±1 Mb) | Enrichment survives all four checks → robust decidual signal (positive result linking to Paper 7) |
| G5 | The cell-context profile is specific to endometriosis | FinnGen leiomyoma, adenomyosis, I9_HYPTENS, T2D | Same pipeline | Profile shared with leiomyoma → specificity withdrawn from G1–G4 |

Exploratory (not in primary families): menstrual/shedding-state enrichment (only 2 menstrual
biopsies in HECA; none in GSE266265) — reported descriptively if any atlas supports it.

## Fixed analysis rules

- LDSC: HapMap3 SNPs, no MHC; 1000G Phase 3 EUR references (GRCh37); FinnGen and the
  harmonised GWAS mapped by rsID. LDSC-SEG uses baseline v1.2 + all-genes control
  annotation; baselineLD v2.2 as sensitivity.
- MAGMA: 1000G EUR reference; window 35 kb up / 10 kb down; sensitivity 0/0 and 100/100.
- Gene sets: top 10% specificity genes per cell type per context; donor-level pseudobulk;
  ≥ 50 cells and ≥ 3 donors per cell type per context; for GSE266265 specificity is computed
  on non-ENDO donors only (all donors as sensitivity) so disease cannot shape the sets.
- Multiple testing: Holm within each hypothesis' primary family; everything else FDR and
  labelled exploratory.
- Negative controls: matched random gene sets (length, expression); I9_HYPTENS, T2D.
- Pivot rule (pre-registered): if no endometrial or lesion context is enriched in G1/G4,
  reframe as "endometriosis genetic risk acts outside the endometrium".

## Compute

MacBook Air M-series, 24 GB RAM, ~320 GB free. Separate Python env for LDSC
(numpy 1.23 / pandas 1.5); main env for everything else. No full atlas loaded in RAM.
