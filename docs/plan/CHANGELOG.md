# Plan change log

- 2026-11-05 v0.1 draft created (pre-preflight).
- 2026-11-05 v0.2 frozen after preflight v001 (before any analysis):
  - FinnGen release fixed to R13 (manifest reachable; R12 not needed).
  - Subtype list fixed with case counts; intestine, fallopian tube, skin scar, NOS removed
    for power/heterogeneity; C3_OVARY_ENDOMETROID removed (ovarian cancer, 208 cases).
  - Sample-overlap rule added: FinnGen is part of GCST90205183, so no replication claims;
    primary GWAS assigned per hypothesis.
  - G4 rewritten: HECA already reported decidual-stroma and macrophage enrichment
    (fGWAS, eutopic only). G4 is now a pre-registered stress test of that result.
    Menstrual-state enrichment moved to exploratory (2 menstrual biopsies in HECA).
  - G2 design fixed to GSE141549 GPL13376, paired lesion vs own eutopic, medication covariate.
  - Non-reproductive control traits fixed: I9_HYPTENS, T2D.
- 2026-11-05 v0.2.1 clarifications fixed before the first analysis run (v003):
  - LDSC implementation: belowlab "2-to-3" Python 3 port of bulik/ldsc (= PyPI ldsc 2.0.1),
    accepted only together with the result of its bundled upstream nose test suite (logged).
  - LDSC inputs: HapMap3 SNPs without MHC (hm3_no_MHC.list); alleles from 1000G Phase 3 EUR
    .bim files; regression LD scores = 1000G Phase 3 LDscore; weights = hm3_no_MHC weights.
  - G3 primary family: all pairs among ENDO_OVARY, ENDO_PERITONEUM, ENDO_DEEP, ENDO_RECTOVAG
    whose observed-scale h2 z >= 4; test = two-sided Wald test of rg = 1; Holm across pairs;
    "rg < 1 supported" requires Holm p < 0.05 and rg < 1. ENDO_STAGE34 is secondary.
  - Heritability scale: observed scale for all traits; liability scale only for all
    endometriosis (population prevalence 0.10), descriptive.
  - G5 at this stage is descriptive (rg with adenomyosis, leiomyoma, hypertension, T2D);
    the G5 specificity test remains the cell-context comparison.
  - GCST90205183 munged with total N = 474,160 (no per-SNP N); used for rg/enrichment, not
    for liability-scale h2.
- 2026-11-05 v004 (technical, no rule change): v003 stopped at munge_sumstats (Python 3 port
  reads a gzipped header as bytes). v004 passes plain-text input; analysis rules unchanged.
  MAGMA download URLs on ctg.cncr.nl return an HTML page; v004 discovers links from it.
- 2026-11-05 v0.2.2 (after G3 part A result, before any further analysis):
  - X2 registered as exploratory: case-case GWAS ovarian vs peritoneal (CC-GWAS); not a primary
    test and cannot rescue G3 part A.
  - G3 part A recorded as NOT SUPPORTED in docs/results/RESULTS_LOG.md.
  - v005 (technical): MAGMA from surfsara share links; MAGMA gene analysis of GCST90205183 used
    only as a sanity check against known loci (no new claim).
- 2026-11-05 v0.2.3 gene-set construction rules (fixed before building any set; v006):
  - Source: GSE266265 raw UMI counts (Paper 8 v2.32c export, read-only); doublets
    (scDblFinder.class) and unannotated cells removed; pseudobulk = sum of UMIs per sample x
    AnnotationRefined.
  - Only non-endometriosis donors define specificity (disease cannot shape the sets).
  - A donor-level group counts if >= 10 cells; a cell type is eligible with >= 3 donors and
    >= 50 cells in total.
  - Expression = CPM per donor group, averaged over donors with equal weight; universe = genes
    mapped to a unique autosomal NCBI37.3 (MAGMA) gene with mean CPM >= 1 in at least one type.
  - Specificity S = mean CPM in type / sum over types; set = top 10% of the universe by S.
  - Family A: AnnotationRefined types (all phases pooled). Family B: AnnotationMain x
    MenstrualCyclePhase (5 levels). ALL_EXPRESSED = control set for LDSC-SEG.
  - MAGMA built from source on arm64 with Eigen vectorization disabled (numerically identical
    algorithm, scalar code path).
- 2026-11-05 v0.2.4 (after v006 data inspection, before any enrichment analysis):
  - Control label in GSE266265 is "CTL" (statuses: Endo / CTL); v006 produced no sets because
    the label filter did not match. Fixed; rules otherwise unchanged.
  - Family B changed to AnnotationMain x MenstrualCyclePhase_main (3 levels): CTL donors have only
    2 donors each in secretory_early and secretory_mid, below the >= 3 donor rule. The 5-level
    version with all donors is kept as an exploratory sensitivity family (SB). A refined-type
    family with all donors (SA) is also built as sensitivity (gene sets only for now).
  - The "counts" slot of GSE266265 is not integer-valued (max 849; likely ambient-corrected).
    CPM-based specificity is unaffected; noted as a limitation.
  - GSE266265 labels follow Tan et al. 2022 nomenclature (eF1-4, dS1-4, fib C7, Prv-*, MPhi1-5,
    mesothelial), transferred by the original authors with Symphony.
  - LDSC-SEG settings: gene window +-100 kb, LD window 1 cM, 1000G Phase 3 EUR, HapMap3 print SNPs,
    baseline v1.2 + family-specific ALL_EXPRESSED control, hm3 no-MHC weights; one-sided
    coefficient p-values; BH within family and trait (descriptive).
  - G4 primary test (check 1 of 4): GCST90205183, family A_refined_CTL, decidual stromal types
    {dS1-myofibroblast, dS2, dS3, dS4}; Holm over the four. "Decidual enrichment replicated" if any
    Holm p < 0.05 (then G4's bold hypothesis is not supported at this check); otherwise check 1
    supports G4. Secondary: stromal@secretory in family B.
  - Traits run through LDSC-SEG: GCST90205183 (primary), FinnGen ENDO_ALL, LEIOMYOMA, ADENOMYOSIS,
    HYPERTENSION, T2D (G5), ENDO_OVARY, ENDO_PERITONEUM, ENDO_DEEP (G3 part B, descriptive here).
- 2026-11-05 v008 (technical, no rule change): v007 stopped before any enrichment result
  (ldsc --l2 could not parse gzipped .annot files; MAGMA linked with older SDKs but the version
  check rejected it). v008 writes plain-text annotations (full, thin as fallback), probes chr22
  first, and keeps the first MAGMA binary that runs. Gene sets and G4 rule unchanged (v0.2.4).
- 2026-11-05 v009 (technical, no rule change): ldsc --l2 --annot requires the full annotation
  file name (v007/v008 passed a prefix). MAGMA v1.10 linked with MacOSX26.sdk; GCST90205183
  gene analysis started in the background as the pre-specified sanity check.
- 2026-11-05 v010 (technical, no rule change): ldsc --print-snps fails in the Python 3 port
  (df.iloc with a boolean Series). LD scores are computed for all 1000G SNPs and restricted to
  HapMap3 rows when per-set files are written; M / M_5_50 are computed over all reference SNPs
  in both cases, so the regression inputs are equivalent.
- 2026-11-06 v011 (technical, no rule change): v010 LD-score workers ran ~17 h (BLAS thread
  oversubscription: 4 processes x all cores). v011 limits each worker to 2 BLAS threads.
  The Python 3 LDSC port writes .l2.ldscore.gz correctly but crashes writing .l2.M / .l2.M_5_50
  (empty M file, rc=1). M = number of SNPs per annotation; M_5_50 = same over SNPs with MAF > 0.05
  (1000G Phase 3 EUR .frq files, the reference genotypes LDSC uses). v011 rebuilds both and checks
  the MAF>0.05 SNP count against the distributed baseline v1.2 M_5_50 for chr22.
- 2026-11-06 v012 (technical, no rule change): low-heat mode for a fanless laptop (1 LDSC worker,
  2 BLAS threads, nice 10; MAGMA sanity run after LDSC-SEG). Results are identical; only runtime changes.
- 2026-11-06 v0.2.5 (after G4 check 1, before v013):
  - G4 check 1 recorded (RESULTS_LOG): NOT REPLICATED at Holm < 0.05 (dS4 Holm 0.061, dS3 0.085).
  - G4 check 2 rule (fixed now): MAGMA v1.10 gene-property analysis of GCST90205183 gene results
    (window 35/10 kb, 1000G EUR); covariates = continuous specificity of each A_refined_CTL type
    (S = mean CPM in type / sum over types) with Average = log2(mean CPM + 1) as conditioning
    covariate (condition-hide=Average, direction=greater); Holm across dS1-myofibroblast, dS2,
    dS3, dS4; "replicated" if any Holm p < 0.05. Family B reported descriptively.
  - Method positive control (descriptive): hypertension should rank vascular smooth-muscle /
    perivascular / arterial types at the top (as seen with LDSC-SEG).
  - Exploratory hypothesis H-sub (generated from v012 LDSC-SEG, NOT a primary test): ovarian
    endometriosis loads on decidual/secretory stroma (dS4) while peritoneal and deep disease load on
    lymphatic/vascular endothelium and mesothelium. Any support from MAGMA on the same FinnGen data
    is not independent confirmation.
- 2026-11-06 v0.2.6 (after G4 check 2, before v014):
  - G4 check 2 recorded: NOT REPLICATED (MAGMA gene-property; best dS3 Holm 0.20) while the
    hypertension positive control passed (5 vascular types q < 0.05).
  - H-sub recorded as not consistent across methods (deep disease: mesothelial/LEC in LDSC-SEG,
    stromal fibroblasts in MAGMA).
  - G4 check 3 rule (fixed now): HECA celltype labels (Marečková 2024 atlas), donor = genotype,
    donors with Endometrial_pathology == C only, Binary Stage != Hormones, superficial and
    whole-uterus biopsies included; same pseudobulk/specificity thresholds as v0.2.3; primary
    method = MAGMA gene-property (passing positive control); GCST90205183; Holm over
    {dStromal_early, dStromal_mid, dStromal_late}; "reproduced" if any Holm p < 0.05.
    LDSC-SEG on HECA sets is secondary (later run).
  - Positive controls (descriptive): hypertension -> Arterial/Venous/perivascular/uSMCs in top 5;
    leiomyoma -> uSMCs (myometrium present only in whole-uterus biopsies).
- 2026-11-06 v0.2.7 (after G4 check 3, before any G1 data processing or genetics):
  - G4 check 3 recorded: NOT REPRODUCED (HECA dStromal_mid rank 1/32, Holm 0.245) with both positive
    controls passing on HECA (leiomyoma -> uSMCs rank 1, q = 0.001; hypertension -> 5/5 vascular).
    Check 4 (leave-one-locus-out) is moot for significance; it may be run descriptively on the
    decidual tilt. Pending: confirm which GWAS / annotation HECA's fGWAS used (v015 retrieves
    their repository) before wording the claim.
  - X3 registered (exploratory): ovarian-subtype tilt towards perivascular / smooth-muscle programs
    (both atlases); to be tested conditional on the leiomyoma genetic component.
  - G1 processing rules (GSE213216): exclude Visium/non-scRNA libraries; cells with >= 200 genes and
    < 25% mitochondrial counts; samples with >= 200 QC cells; per-sample clustering (2,000 HVG,
    30 PCs, Leiden resolution 1.0; k-means fallback); cluster labels by Spearman correlation of
    cluster log2(CPM+1) with GSE266265 Tan-nomenclature centroids (top 30 specific genes per type,
    all donors); independent canonical-marker lineage reported as agreement; pseudobulk = sum of raw
    counts per sample x label.
  - G1 primary test (fixed now): GCST90205183; contrast peritoneal lesion vs eutopic endometrium within
    each of 5 lineages (stromal, epithelial, endothelial, myeloid, lymphocyte; lineage = reference
    label's AnnotationMain); per lineage, patient-level log2(CPM+1) pseudobulk (>= 10 cells), >= 3
    patients per context; covariate = mean lesion - mean eutopic log2(CPM+1); conditioning covariate =
    mean log2(CPM+1) over both contexts; MAGMA gene-property direction=greater; Holm over lineages.
    "G1 supported" if any lineage has Holm p < 0.05 with beta > 0. The same model with
    direction=less is the pre-specified alternative (source-defect model). Endometrioma vs eutopic and
    FinnGen subtypes are secondary; mesothelial (lesion-only) cells are descriptive.
- 2026-11-06 v0.2.8 (before v016):
  - HECA's published fGWAS used the same GWAS (GCST90205183, GRCh38) with per-cell-type mean
    expression at gene TSS (their repository, notebooks Endometriosis_fGWAS1-4). Our G4 claim is
    therefore method-level: same GWAS, same atlas labels, decidual enrichment not reproduced by
    LDSC-SEG / MAGMA gene-property while positive controls pass.
  - G4 check 4 (fixed now; mechanism for the discrepancy): (a) mean specificity percentile of the
    Bonferroni-significant MAGMA genes of GCST90205183 in each HECA type, empirical p from 10,000
    random gene sets of equal size (seed 20261106); (b) MAGMA gene-property on HECA types after
    removing all genes within +-1 Mb of those genes. Prediction if a few loci drive fGWAS: dStromal
    high in (a), dStromal beta shrinks / rank drops in (b). Reported descriptively (no new primary
    claim); the result informs the interpretation of the HECA-vs-polygenic discrepancy.
  - v016 (technical): exact kNN (sklearn) instead of pynndescent (numba compile failure for samples
    > 4,096 cells); per-sample cache; non-tar.gz archive members extracted. Processing rules unchanged.
- 2026-11-06 v0.2.9 (after GSE213216 preparation QC, before any G1 genetics):
  - G4 check 4 recorded (RESULTS_LOG): top loci are not dStromal-specific; removing +-1 Mb around
    them barely changes the dStromal tilt -> "a few strong loci" does not explain HECA vs polygenic.
  - Label transfer quality: lineage agreement 90.6% (cell-weighted), but fine-label margins are tiny
    (median 0.012), so only lineage-level analyses are interpretable (as already specified for G1).
  - G1 filter (new, data-quality motivated, decided before genetics): only clusters whose reference
    lineage and canonical-marker lineage agree enter G1 (removes 9.4% of cells, mainly reference
    "lymphocyte" clusters with stromal marker profiles).
  - G1 details fixed: patient-level profiles sum all agreeing clusters of a lineage over that
    patient's samples of the same context (>= 10 cells); genes = unique autosomal NCBI37.3 genes with
    CPM >= 1 in either context; conditioning covariate AVG = mean of the two context means.
  - Secondary G1 analyses (descriptive): matched subtypes (peritoneal lesion contrast x ENDO_PERITONEUM;
    endometrioma contrast x ENDO_OVARY), ENDO_ALL, specificity controls (HYPERTENSION, LEIOMYOMA),
    endometrioma vs eutopic (GCST), paired patients (both contexts) for GCST.
  - G4 sensitivity: HECA all donors (any pathology, hormones included), same specificity rules.
- 2026-11-06 v018 (technical, no rule change): MAGMA v1.10 accepts only direction=greater (and
  two-sided); the pre-specified one-sided "less" test is implemented as direction=greater on DN = -D
  (identical hypothesis). v017 produced no G1 statistics because the "less" call failed.
- 2026-11-06 v0.2.10 (after the G1 primary run v018; before v019; no GTEx or gnomAD data touched yet):
  - G1 recorded (RESULTS_LOG): ALTERNATIVE by the pre-registered rule (lymphocyte, Holm p_less 0.017);
    no lesion-upregulated enrichment in any lineage. Post hoc observation: the HYPERTENSION control
    shows the same eutopic-up direction in 4/5 lineages (Holm <= 0.03), with effects as large as for
    endometriosis; LEIOMYOMA does not (registered G5 falsifier "profile shared with leiomyoma" not
    met); the subtype-matched ENDO_PERITONEUM shows nothing. The G1 alternative is therefore
    reported as not endometriosis-specific. Analyses G1-S below are post hoc diagnostics
    (descriptive, labelled as such), not new confirmatory tests.
  - G1-S diagnostics (post hoc), on the v018 covariate files (data/processed/magma_covar/g1):
    S1: second unrelated control trait T2D (FinnGen R13; MAGMA gene analysis as in v013).
    S2: generic gene-property explanation: conditioning set extended to AVG + gnomAD v2.1.1 LOEUF +
        GTEx v8 expression breadth (tau over the 52 G6 tissues); genes with all three values; the
        AVG-only model is refitted on the same genes for comparison. Interpretation fixed now:
        "generic gene-property explanation supported" if in the lymphocyte contrast both
        GCST90205183 and HYPERTENSION betas shrink by >= 50% and reach p_less >= 0.05.
    S3: GCST90205183 conditional on the HYPERTENSION gene z (and, separately, on the T2D gene z).
    S4: label composition of the lymphocyte lineage per context (descriptive).
  - G6 registered (core; "outside the uterus" arm of the central question):
    GTEx v8 gene median TPM, 52 tissues (the two "Cells -" lines excluded); covariate per tissue =
    log2(TPM + 1); conditioning covariate Average = mean over the 52 tissues (condition-hide);
    MAGMA gene-property, direction=greater; genes = unique autosomal NCBI37.3 symbols with
    TPM >= 1 in >= 1 tissue. Primary trait GCST90205183. Reproductive set R = {Uterus, Ovary,
    Vagina, Fallopian Tube, Cervix - Ectocervix, Cervix - Endocervix}. Verdict: "REPRODUCTIVE" if any
    R tissue has Holm (over the 6 R tissues) p < 0.05; otherwise "NON-REPRODUCTIVE" if any other
    tissue has BH q < 0.05 over all 52; otherwise "NO TISSUE CONCENTRATION".
    Secondary (descriptive): ENDO_ALL, ENDO_OVARY, ENDO_PERITONEUM, ENDO_DEEP, LEIOMYOMA,
    HYPERTENSION, T2D. Positive controls: HYPERTENSION -> an Artery tissue in the top 3;
    LEIOMYOMA -> Uterus in the top 3.
    Context: Rahmioglu 2023 reported endometrium / uterus / smooth-muscle enrichment for genes near
    the 49 index SNPs (TissueEnrich); G6 asks whether the genome-wide polygenic signal follows
    (parallel to G4 for cell types). Limitation fixed now: GTEx "Uterus" is mostly myometrium.
- 2026-11-06 v0.2.11 (after an external review of the plan; before the v019 log was seen; no new
  genetic results seen):
  - Erratum: GSE213216 is Fonseca et al. 2023, Nat Genet 55:255-267 (Lawrenson lab; >370,000
    cells; 17 patients with endometriosis + 4 unaffected), not "Tan 2022" as written in PLAN_v0.2
    and in the v015-v018 script headers. Tan et al. 2022 (Nat Cell Biol) is GSE179640; only the
    label nomenclature of our reference (GSE266265 annotations from Paper 8) follows Tan.
    No analysis is affected.
  - Dataset overlap (documented): HECA integrates Wang 2020 (GSE111976), Garcia-Alonso 2021
    (E-MTAB-10287), Tan 2022 (GSE179640), Lai 2022 (GSE183837), Fonseca 2023 (GSE213216),
    Huang 2023 (GSE214411) and new data. GSE266265 (60 donors, GEO 2026) is not part of HECA.
    G1 (GSE213216) and G4 (HECA, GSE266265) are not presented as independent confirmations of each
    other; all analyses use one GWAS (GCST90205183, which contains FinnGen), so no genetic
    replication is claimed anywhere.
  - Literature overlap (docs/plan/LITERATURE_OVERLAP.md): Koller et al. 2026 Nat Genet
    (multi-ancestry, 105,869 cases, 80 loci; TWAS in 7 GTEx tissues; colocalization; includes
    23andMe) and Guare et al. 2025 (GBMI; scDRS on HECA, endometriosis vs non-endometriosis donors).
    Consequences fixed now: G6 is re-classified as exploratory (tissue-level context only);
    G2 keeps only its distinct question (sign concordance of genetically predicted expression with
    lesion-vs-eutopic vs patient-vs-control expression) and will use published TWAS results where
    available; the paper's core is a pre-registered stress test of cell-context localization
    (G4, G1 with trait-specificity controls, explicit detection limits). Whether the published
    Koller version contains HECA S-LDSC is unverified (absent from its medRxiv version) and must be
    checked from the full text before G4 is worded.
  - Interpretation rule (fixed now): non-significant enrichment is reported as "not detected at the
    power of this GWAS / method", with detection limits, never as "risk does not act there".
  - G4 check 5 (annotation design as the mechanism of the HECA-vs-polygenic discrepancy), on the
    v014 HECA control-donor pseudobulk, MAGMA gene-property per type:
    M1 specificity | Average (v014 primary, refit); M2 log2(mean CPM + 1) | Average;
    M3 log2(mean CPM + 1) without conditioning (closest analogue of the fGWAS mean-expression
    annotation). Traits GCST90205183, HYPERTENSION, LEIOMYOMA (T2D if available).
    Reading fixed now: "annotation design explains the discrepancy" if a dStromal type reaches
    Holm p < 0.05 (over the 3 dStromal types) under M3 but not under M1, and the M3 top-3 types of
    GCST90205183 and HYPERTENSION share >= 2 types (non-specific). Descriptive.
  - Detection limits (fixed now): for G4 (HECA M1, 3 dStromal types) and G1 (5 lineages, both
    directions): minimal detectable effect at 80% power, one-sided alpha = 0.05/k (k = 3 or 5),
    MDE = (z_(1-alpha) + z_0.8) x SE, in BETA and BETA_STD units, next to the observed estimate
    (95% CI) and the positive-control effects (LEIOMYOMA -> uSMCs; HYPERTENSION -> best
    vascular type) in their own GWAS.
  - Fonseca annotation check (data quality, before any re-analysis): the authors' GEO-to-patient /
    class table (github.com/lawrenson-lab/AtlasEndometriosis) is compared with our GEO-derived
    patient and context assignment; the repository file list is inspected for per-cell
    annotations. Planned G1 sensitivity (rule fixed now): the G1 primary model re-run with the
    authors' cell annotations (their lineage labels replace label transfer; same thresholds),
    once their Seurat objects are available.
- 2026-11-06 v0.2.12 (after the second part of the external review; before v019/v020 logs were seen):
  - Ochoa et al. 2025 (Adv Sci; bioRxiv 2024.11.21.624726) added to LITERATURE_OVERLAP: scDRS on
    GSE213216 lesion cells only (200,798 cells), GWAS 23,492 cases; M2 macrophages p = 0.013,
    endothelial 0.044, dendritic 0.050 (uncorrected Monte Carlo p); IL1 follow-up experiments.
    "First link of endometriosis genetics to lesion cells" is therefore not a claim of Paper 9.
  - Koller et al. 2026 resolved from the author manuscript (PMC, 2026-05-22; Nat Genet 58:1051-1061):
    S-LDSC on HECA cell-type marker annotations (|log2FC| > 1.2, FDR < 0.01, cases and controls
    separately) gave only "nominally significant enrichment for decidual stromal cells and SOX9
    functionalis" (their Supplementary Tables 9-10) with 105,869 cases. Their summary statistics
    are on Zenodo (doi 10.5281/zenodo.18983492); 23andMe statistics are restricted. The v021 audit
    lists the Zenodo record (A0) before any decision to use it.
  - Clarification of G1 design: the G1 statistic is already a direct contrast (covariate
    D = lesion - eutopic per gene, one test per lineage), not two separate p-values. A direct
    difference test is still required for cross-trait comparisons (endometriosis vs control
    trait); it is added below with the z-correlation taken from the LDSC cross-trait intercept.
  - Wording rule for G4 (fixed now): "statistically significant prioritization of decidual stroma
    was not reproduced under the pre-specified criteria, although individual decidual states rank
    high; ranks are not evidence of an effect" (no "refuted").
  - Audit before any new hypothesis (v021, descriptive; results decide whether earlier outputs
    are kept):
    A1 sumstats: real rows / non-missing rows per munged file (munge writes every merge-allele SNP,
       so "Writing summary statistics for N SNPs" counts NA rows; the v004 log reported that N);
    A2 alleles / build: sign concordance of genome-wide significant GCST90205183 SNPs with FinnGen
       ENDO_ALL; genome-wide Z correlation; MAGMA SNP mapping counts per trait (rsID-based, so the
       GRCh37 / GRCh38 difference between the original and harmonised GCST files does not enter);
    A3 LDSC h2 diagnostics per trait (intercept, ratio, mean chi2, lambda GC, N), flag ratio > 0.3;
    A4 provenance of the non-integer GSE266265 counts (fraction integer, fractional parts,
       per-cell sums vs nCount_RNA if exported);
    A5 GSE266265 leave-one-donor-out: MAGMA gene-property for the dS / stromal types of families
       A_refined_CTL and B_main_phase3_CTL with each control donor removed in turn;
    A6 G1 sensitivity: (i) all clusters (no lineage-agreement filter), (ii) >= 50 cells per patient
       profile, (iii) leave-one-patient-out (lymphocyte; GCST90205183 and HYPERTENSION),
       (iv) patient-label permutation null for the lymphocyte contrast (100 permutations, both
       traits; calibration of the parametric MAGMA p-value when D comes from 5-12 patients);
    A7 direct cross-trait difference of the G1 betas (GCST90205183 vs HYPERTENSION, vs T2D if
       available): z = (b1 - b2) / sqrt(se1^2 + se2^2 - 2 rho se1 se2), rho = LDSC gcov_int.
  - Framing: PLAN_v0.3_DRAFT.md (not frozen) proposes the question "How robust is cell-context
    prioritization of inherited endometriosis risk across GWAS, single-cell references and
    analytical frameworks?"; it is frozen only after the v021 audit and with the original
    hypotheses (v0.2) kept and labelled separately from post-review ones.
- 2026-11-06 v0.2.13 (after the v019 log; before v021 / v020 logs; no X4 data touched):
  - G6 and G1-S recorded (RESULTS_LOG). G6 positive control for hypertension failed; leiomyoma
    passed; G6 stays exploratory.
  - X4 registered (exploratory, post hoc: formed after seeing G6 strong reproductive-tissue and
    G4 weak cell-type results). "Organ-identity hypothesis": endometriosis risk genes are
    reproductive-tract (Muellerian) identity genes expressed across many endometrial cell types,
    which would explain strong tissue-level but weak cell-type-level enrichment.
    Definitions (GTEx v8, the G6 matrix): MULL = mean log2(TPM+1) over {Uterus, Fallopian Tube,
    Cervix - Ectocervix, Cervix - Endocervix, Vagina}; SMOOTH = mean over {Colon - Sigmoid,
    Esophagus - Muscularis, Bladder, Artery - Tibial}; Average = mean over the 52 tissues.
    HECA within-endometrium breadth: tau over the HECA control-donor types (log2 mean CPM + 1);
    genes in both universes split into tau tertiles (T1 broad ... T3 cell-type specific).
    X4a: MULL | Average + SMOOTH (and Uterus alone | Average + SMOOTH), direction=greater.
    X4b: X4a model within each tau tertile (disjoint gene subsets); BETA_STD compared T1 vs T3.
    X4c: HECA dStromal_mid / dStromal_early specificity | Average_HECA, then additionally | MULL.
    Traits: GCST90205183 (primary), ENDO_ALL, ENDO_OVARY, ENDO_PERITONEUM, ENDO_DEEP, LEIOMYOMA,
    HYPERTENSION, T2D.
    Reading fixed now: "Muellerian-specific, not smooth muscle" if GCST90205183 X4a MULL p < 0.05;
    "endometriosis-specific" additionally requires HYPERTENSION and T2D X4a p >= 0.05;
    "organ-identity pattern" if GCST90205183 X4b T1 p < 0.05 and BETA_STD(T1) > BETA_STD(T3) with
    one-sided p < 0.05 for z = (b1 - b3) / sqrt(se1^2 + se3^2). Any support is exploratory and
    needs a second GWAS (e.g., Koller 2026 Zenodo statistics, if suitable after A0) to count.
- 2026-11-06 v0.2.14 (after v020, v021, v022 logs; before any Koller 2026 data is read):
  - v020-v022 recorded (RESULTS_LOG). G1 statistics: parametric MAGMA p-values for contrasts
    estimated from few donors are anti-conservative (permutation t SD 1.37); from now on every
    donor-level contrast is also reported with a patient-label permutation p (B >= 100).
  - Second GWAS (larger, partially overlapping; not a replication): Koller et al. 2026
    endometriosis_EUR_wo_23andMe (Zenodo 10.5281/zenodo.18983492, CC-BY-NC-4.0). Shares UKB,
    FinnGen, EstBB and IEGC participants with GCST90205183, so agreement shows robustness to GWAS
    choice, not independence. Processing: rsID-based (chr:pos mapped to rsID via the GCST
    harmonised GRCh38 file or the 1000G GRCh37 .bim, whichever matches); MAGMA gene analysis as
    v013 (35/10 kb, g1000_eur); munge as v004; LDSC rg with GCST90205183 and ENDO_ALL (overlap).
    Re-run unchanged on it, rules fixed now: G6 (GTEx, 52 tissues), X4a, X4b, G4 HECA M1 and M3,
    G1 primary (lesion vs eutopic, 5 lineages) with permutation p for lymphocytes.
    Reading: a finding "holds with the larger GWAS" if same sign and p < 0.05 (Holm where the
    original rule used Holm); G4 counts as reproduced only if a dStromal type reaches Holm < 0.05
    under M1. Any non-EUR file (EAS, AFR, AMR) is not used yet (needs ancestry-matched LD).
  - .gitignore: references/** and *.pdf, *.xlsx, *.docx are never committed (publisher files).
- 2026-11-06 v0.2.15 (after the v023 log; before v024):
  - v023 recorded. G4 status changes: with the larger (partially overlapping) GWAS the
    pre-specified G4 criterion is met (dStromal_mid Holm 0.037). Specificity is the open question.
  - G4-S (post hoc, descriptive; fixed now, HECA control donors, M1 = specificity | Average), on
    the genes that have all covariates, traits KOLLER_EUR and GCST90205183:
    S5a dStromal_mid (and early / late) | Average + gene z of T2D, of HYPERTENSION, of LEIOMYOMA
        (one at a time). Reading: "not explained by generic gene-level signal" if dStromal_mid stays
        p < 0.05 under all three.
    S5b dStromal_mid | Average + specificity of ePV_2, Arterial and uSMCs. Reading: "separable from
        perivascular / smooth-muscle programs" if p < 0.05.
    S5c dStromal_mid | Average + MULL (GTEx Muellerian expression).
    S5d M1 panel for every available trait (dStromal p, Holm over 3, rank among 32).
  - G1-S (post hoc): patient-label permutation p (B = 100) for the epithelial contrast, traits
    KOLLER_EUR, GCST90205183, T2D.
- 2026-11-06 v0.2.16 (after the v024 log; before v025):
  - v024 recorded. G4 is recorded as "criterion met with the larger GWAS but not endometriosis-
    specific (T2D meets it; not separable from perivascular / SMC programs)".
  - Option A component "Ochoa 2025 claim" (fixed now): scDRS on GSE213216 cells.
    Cells: the v016 per-sample pipeline is re-run to give per-cell labels (cluster label
    transfer as before; lineage-agreement filter); contexts: lesion = peritoneal lesion +
    endometrioma (as Ochoa), eutopic separately; cells subsampled to <= 150 per sample x label
    (seed 20261106) to fit memory. Expression: counts normalized to 1e4, log1p; scDRS preprocess
    with covariates const + n_genes; gene sets: top 1,000 MAGMA genes by z (z capped at 10 as
    weights), traits GCST90205183, KOLLER_EUR, HYPERTENSION, T2D, LEIOMYOMA; n_ctrl = 1,000,
    seed 0; group test = scDRS association test (0.95 quantile of normalized scores in the
    group vs control gene sets, Monte Carlo p).
    Groups: lineage (5) and fine labels with >= 100 cells. Macrophage labels = labels starting
    with "MPhi" / "M$\Phi$" in the Tan nomenclature (Ochoa's "M2 macrophages" have no exact label
    here).
    Reading: "Ochoa claim reproduced" if the myeloid lineage or any macrophage label has
    Holm (over macrophage labels + myeloid lineage) MC p < 0.05 in lesion cells for GCST90205183
    or KOLLER_EUR; "endometriosis-specific" additionally requires HYPERTENSION and T2D MC p >= 0.05
    for the same group. Eutopic results are descriptive (new; not tested by Ochoa).
- 2026-11-06 v026 (technical, no rule change): the v025 scdrs install pulled scikit-misc (cp39 wheel
  built against numpy 1.x) into the main .venv, which runs numpy 2.x -> import error. scDRS now runs
  in a separate .venv_scdrs with a pinned numpy-1.26 stack; the main .venv is only probed and
  repaired if pipeline modules fail. The v025 analysis (rules v0.2.16) is re-run unchanged.
  Per-cell labels are recomputed with scanpy 1.10.3 in that environment (same processing rules).
- 2026-11-06 v0.2.17 (after the second external review and the choice "A main + B exploratory";
  before v027; v025/v026 scDRS results not yet seen):
  - PLAN_v0.3.md frozen (replaces PLAN_v0.3_DRAFT.md).
  - A7 (v021) cross-trait difference tests are descriptive only: gcov_int is not the error
    correlation of MAGMA gene-property coefficients. Replaced by delete-one-block jackknife (A2).
  - v027 analyses (fixed now; HECA control donors, M1 = specificity | Average):
    R1 correlation (Pearson, Spearman) among specificity vectors of dStromal_early/mid/late,
       ePV_1a, ePV_1b, ePV_2, mPV, Arterial, Venous, uSMCs (those present); VIF of dStromal_mid
       given ePV_2 + Arterial + uSMCs and given all other types.
    R2 pairwise conditional matrix: p of type i | Average + type j for the focus set
       {dStromal_mid, ePV_2, ePV_1b, Arterial, uSMCs, mPV}; traits KOLLER_EUR, GCST90205183,
       HYPERTENSION, T2D, LEIOMYOMA.
    R3 delete-one-block jackknife, 50 contiguous gene blocks (genome order): betas of dStromal_mid,
       ePV_2, Arterial, uSMCs (marginal | Average) and of the G1 contrasts (lymphocyte,
       epithelial; D | AVG) for KOLLER_EUR, GCST90205183, HYPERTENSION, T2D; jackknife SEs for
       within-trait differences (dStromal_mid - ePV_2, standardized betas) and cross-trait
       differences (block-paired).
    R4 simulations with known truth: null gene z = KOLLER_EUR gene z circularly shifted along the
       genome (random offset; keeps local correlation, breaks alignment); signal added as
       b * standardized specificity of the true type; scenarios null, dStromal_mid, ePV_2,
       Arterial, uSMCs at 1x and 2x the observed KOLLER dStromal_mid standardized effect;
       40 replicates each; read-outs: P(dStromal_mid top), P(dStromal_mid Holm-3 < 0.05),
       P(true type top), P(true type p < 0.05), false-positive rate under the null.
    R5 calibration under a true null: random splits of patients within one context (eutopic
       stromal, eutopic epithelial, lesion stromal), D | AVG, KOLLER_EUR, 200 splits each;
       read-out: fraction of parametric p < 0.05 (nominal 5%).
    B0 discovery only: list CELLxGENE collections / datasets matching "fallopian" (no download).
    K0 if Koller supplementary tables exist in references/, print Supplementary Tables 9 and 10.
- 2026-11-06 v0.2.18 (after the v025 log; before v028):
  - v025 recorded (A4: Ochoa claim not reproduced on our labels; hypertension shows the myeloid
    association instead).
  - v028 Ochoa-faithful scDRS (fixed now; .venv_scdrs; otherwise v025 rules):
    contexts: lesion_ochoa = peritoneal lesion + the 4 "no endometriosis detected" peritoneum
    samples (re-classified as peritoneal endometriosis by Ochoa) + endometrioma; lesion_ours as in
    v025; eutopic. Cells: <= 400 per sample x label (seed 20261106).
    Groups: lineages, labels with >= 100 cells, MACRO_ALL (all MPhi labels), M2like = myeloid cells
    whose mean z-scored expression of CD163, MRC1, LYVE1, F13A1, STAB1, SELENOP, FOLR2 is above
    the myeloid median.
    Tests: scDRS association MC p (1,000 controls) with BH across labels per trait and context;
    donor-level test per group: per sample with >= 10 group cells, z = (mean normalized disease
    score of the group's cells - mean of the same for the 1,000 matched control gene sets) / SD
    of the control means; one-sided one-sample t-test of z > 0 across samples.
    Reading: "Ochoa claim reproduced" if in lesion_ochoa MACRO_ALL, M2like, myeloid or a MPhi label
    has Holm (over these groups) MC p < 0.05 for GCST90205183 or KOLLER_EUR AND donor-level p < 0.05
    for the same group; "endometriosis-specific" additionally requires HYPERTENSION and T2D Holm
    MC p >= 0.05 for that group.
  - Runner: v028 first completes v027 if needed (ordering of installers).
- 2026-11-06 v0.2.19 (after the v027 log; before any fallopian-tube data are downloaded):
  - v027 recorded (RESULTS_LOG).
  - B1 fixed now (exploratory, PLAN v0.3): reference = CELLxGENE dataset
    0436a180-cb44-47ba-8ffa-807b7a468469 (Weigert et al. 2025, Nat Commun; fallopian tube, pre- and
    post-menopausal). Pre-menopausal donors only (a menopause column if present, else age < 50, else
    all, logged). Cell-type column: the author annotation column with most categories among obs
    columns named like author / annotation / celltype / cluster; fallback cell_type. Raw counts
    (raw/X if present, else X). Pseudobulk donor x type; types with >= 3 donors (>= 10 cells each)
    and >= 50 cells; specificity S = mean CPM / sum over types; Average = log2(mean CPM + 1);
    genes = unique autosomal NCBI37.3 symbols present in both this atlas and HECA.
    Compartments: FT_SEC = sum of S over secretory types (labels matching "secret|SE[0-9]|PAX8");
    FT_CIL = sum over ciliated types ("cilia|CE" as a word); HECA_EPI = sum of HECA control-donor
    S over non-ciliated epithelial types (names matching SOX9 | Lumin | Gland | MUC5B | KRT5);
    HECA_CIL = sum over HECA ciliated types (names matching "iliated"). The assignment of every
    label to a compartment is printed before any genetics is run.
    Tests (MAGMA gene-property, direction=greater): (i) FT types | Average_FT (descriptive top list);
    (ii) FT_SEC | Average_FT + Average_HECA + HECA_EPI and HECA_EPI | the same + FT_SEC; FT_CIL vs
    HECA_CIL likewise; (iii) delete-one-block jackknife (50 blocks) of BETA_STD(FT_SEC) -
    BETA_STD(HECA_EPI) and FT_CIL - HECA_CIL from the marginal models (| both Averages).
    Primary trait KOLLER_EUR; secondary GCST90205183, ENDO_OVARY, ENDO_PERITONEUM; controls
    HYPERTENSION, T2D, LEIOMYOMA.
    Reading: "FT-preferential" for a compartment if, for KOLLER_EUR, (ii) p < 0.05 AND (iii)
    difference > 0 with jackknife p < 0.05, AND neither HYPERTENSION nor T2D has (ii) p < 0.05.
    No claim about lesion origin in any case.
  - Runner: v029 completes v028 (Ochoa-faithful scDRS) first if needed.
- 2026-11-07 v0.2.20 (after the v028 log; before v030):
  - v028 recorded (A4 Ochoa: not reproduced; donor-level macrophage excess shared with hypertension
    and T2D). Technical note: future scDRS runs process each context in a separate process.
  - A4 Koller 2026 method-matched reproduction (fixed now; their Methods: HECA scRNA-seq, cases and
    controls separately, each cell type vs all other cell types, |log2FC| > 1.2 and FDR < 0.01,
    endometrial and immune compartments separately, S-LDSC):
    cells: HECA donors of E-MTAB-14039 if a dataset column identifies them (values matching
    "14039" or "Marec"), otherwise all HECA donors; status control = Endometrial_pathology "C",
    case = categories naming endometriosis ("E..." / "endometriosis"); "Binary Stage" = Hormones
    excluded. Compartments: immune = labels matching
    Immune|^uM|uNK|^T_|^B_|Mast|DC|Mono|Macro|Neutro|Plasma|ILC|Lymphoid|Myeloid (not Lymphatic);
    endometrial = all others; the label assignment is printed before testing.
    Markers: per status x compartment, cell types with >= 30 cells, <= 400 cells per type (seed
    20261106); log-normalized; scanpy Wilcoxon each type vs rest; genes = unique autosomal NCBI37.3
    symbols; sets ABS (|log2FC| > 1.2, adjusted p < 0.01; primary, as written) and UP
    (log2FC > 1.2; sensitivity); sets with < 20 genes dropped; control = all tested genes.
    S-LDSC (LDSC-SEG layout, as v012): +-100 kb windows, 1 cM, baseline v1.2, HapMap3 weights,
    --h2-cts with the all-genes control, one-sided coefficient p, BH within family; family = ABS.
    MAGMA gene-set analysis (competitive, condition-hide=ALLGENES) for ABS and UP as a cheap
    complement. Traits KOLLER_EUR (primary), GCST90205183, HYPERTENSION, T2D, LEIOMYOMA.
    Reading: "Koller nominal result reproduced" if for KOLLER_EUR (S-LDSC, ABS) at least one
    dStromal annotation and at least one SOX9 functionalis annotation have p < 0.05; "survives
    correction" if any of them has q < 0.05; "endometriosis-specific" if the reproduced annotations
    have p >= 0.05 for HYPERTENSION and T2D.
  - Runner: v030 completes B1 (v029) first if needed.
- 2026-11-07 v0.2.21 (after the v030 log; before v031):
  - v030 recorded (A4 Koller: nominal result reproduced, not BH-significant; the BH-significant
    annotations under the same pipeline are endothelial / perivascular).
  - A4 follow-up, Koller pipeline (v031). All on the v030 cells (Mareckova, non-hormone, <= 400
    cells per type, seed 20261106) unless stated; markers re-derived with identical code; if the
    ABS sets equal the v030 sets exactly, the v030 LD scores and S-LDSC results are reused.
    K1 cohort audit: donors and samples (first obs column named sample / sample_id / library) by
       status, with and without hormone-treated samples, Mareckova and all HECA; compared with
       Koller's 22 endometriosis / 13 control samples (printed, no test).
    K2 set composition: for every ABS set, the share of genes with log2FC < -1.2 (DOWN share);
       Jaccard overlap between the focus sets and the top sets.
    K3 selection rule: S-LDSC families UP (log2FC > 1.2) and DOWN (log2FC < -1.2), same layout.
       Prediction H-ABS (stated now): endothelial ABS sets (Venous, Arterial) have DOWN share > 0.5;
       in UP, no endothelial set is in the KOLLER_EUR top 6 and a dStromal set is; in DOWN,
       endothelial sets rank in the top 6. H-ABS "supported" if all three hold.
    K4 conditional resolution (ABS and UP families; KOLLER_EUR and GCST90205183): rows = focus sets
       (labels dStromal* / SOX9_functionalis*), columns = the 6 best non-focus sets by KOLLER_EUR
       p in that family; each ordered pair fitted jointly (--h2-cts line "A,ALLGENES,B": coefficient
       of A given B, all genes and baseline). For a pair with both marginal p < 0.05: "A survives"
       if p(A|B) < 0.05; outcome both / A only / B only / neither. Reported: focus sets that survive
       conditioning on every column set, and column sets that survive every focus set.
    K5 sensitivity (ABS, KOLLER_EUR, GCST90205183, ENDO_ALL): (H) hormone-treated samples included;
       (N) all cells (no 400 cap). Same reading as v0.2.20 per variant.
    K6 comparator informativeness: ENDO_ALL (FinnGen R13, same cohort and LD mismatch as the
       control traits) added to every family. If the Koller-reproduced annotations have p >= 0.05
       for ENDO_ALL as well, absence in HTN/T2D is not read as specificity. S-LDSC positive controls
       reported: HYPERTENSION q < 0.05 for a vascular/perivascular/smooth-muscle set; LEIOMYOMA
       q < 0.05 for uSMCs or mPV.
  - Disk: per-family annotation and full LD-score files are deleted after the HapMap3 split
    (v031 families only).
- 2026-11-07 v0.2.22 (after the v029 log; before v032):
  - v029 recorded: compartment step not run (label patterns matched no author label).
  - Amendment to B1 labels (technical; rules otherwise unchanged): the Weigert author labels use
    the prefixes se_ (secretory epithelium) and ce_ (ciliated epithelium); FT_SEC adds "^se_",
    FT_CIL adds "^ce_". The assignment is checked against the CELLxGENE cell_type column (printed:
    share of cells per standard cell type for each compartment).
    Disclosure: this amendment is written after v029 printed the descriptive per-type top-5 lists
    (no se_/ce_ type in any endometriosis top 5; st3 first for KOLLER_EUR and ENDO_OVARY).
    No compartment-level test (ii)/(iii) has been run with any label set.
  - POST HOC (suggested by the v029 list; hypothesis-generating only): FT_STR = sum of S over
    st* types vs HECA_STR = sum over HECA types matching "Stromal|Fibroblast"; same tests (ii),
    (iii) and the same reading as the v0.2.19 compartments, labelled POST HOC in every output.
  - Completion marker printed only when the tests finished.
- 2026-11-07 technical (v033 installer; no rule change): v031 was interrupted by the user during
  the P_UP LD scores; v033 writes the updated docs (v032 results) and resumes v031 unchanged
  (completed LD-score chromosomes and cached S-LDSC results are reused).
- 2026-11-07 v0.2.23 (after the v031 log; before v034):
  - v031 recorded. H-ABS NOT SUPPORTED by its rule (c1 failed on Arterial/ENDO).
  - K7 sign decomposition (new; stated after seeing the v031 marginal UP and DOWN families, so the
    marginal pattern is known; the joint fits are not): for every cell type and status with UP and
    DOWN sets, joint S-LDSC "UP,ALLGENES,DOWN" and "DOWN,ALLGENES,UP" on the v031 LD scores. For
    types with ABS p < 0.05: "UP-carried" if p(UP | DOWN) < 0.05 and p(DOWN | UP) >= 0.05;
    "DOWN-carried" the reverse; "both"; "neither". Traits KOLLER_EUR (primary), GCST90205183,
    ENDO_ALL. H-SIGN "supported" for KOLLER_EUR if every ABS-nominal Venous/Arterial type is
    DOWN-carried and every ABS-nominal dStromal type is UP-carried.
  - K8 (prediction from G4: decidual stroma not separable from perivascular ePV_2): UP sets,
    rows = dStromal sets, columns = Arterial, Venous, ePV_1a, ePV_1b, ePV_2, mPV, uSMCs (both
    statuses where present); both directions. "G4 prediction confirmed" for a trait if
    p(dStromal_mid | ePV_2) >= 0.05 for at least one dStromal_mid x ePV_2 pair; "not confirmed"
    if all such pairs give p < 0.05.
  - No new LD scores; uses the v031 P families.
- 2026-11-07 v0.2.24 (after the v034 log; before v035):
  - v034 recorded. H-SIGN NOT SUPPORTED by its rule (Arterial/ENDO UP-carried); K8: G4 prediction not
    confirmed in S-LDSC up-only sets.
  - S-ABS simulation (v035), known source, MAGMA competitive gene-set analysis (condition-hide
    ALLGENES) of the v031 P-cell ABS / UP / DOWN sets in one model run per replicate.
    Background: KOLLER_EUR gene z circularly shifted (offset uniform in [500, n-500]), seed 20261107.
    Planted truth: standardized HECA (heca_v014) specificity profile of the target type added to z
    at 1x and 2x the R4 effect (0.947 x sd of the dStromal_mid profile). Scenarios: null;
    dStromal_mid x1, x2; Venous x1, x2; ePV_2 x2. 60 replicates each. Observed KOLLER_EUR run once.
    Groups of the top set: dStromal / endothelial (Arterial, Venous) / perivascular+SM (ePV_*,
    mPV, uSMCs) / other. BH within family.
    Readings (primary strength 2x): calibrated if null any q < 0.05 <= 0.10 per family;
    "ABS relabels a decidual source" if, under dStromal_mid truth, P(top ABS set endothelial or
    perivascular/SM) > P(top ABS set dStromal) AND exceeds its own null share by >= 0.10;
    "UP recovers the source" if P(top UP set dStromal) >= 0.5. Venous and ePV_2 truths reported
    descriptively (does each rule recover a true vascular source).
- 2026-11-07 v0.2.25 (after the v035 log; before v036):
  - v035 recorded: null not calibrated (any q < 0.05: ABS 0.17, UP 0.27); ABS relabel NO; UP
    recovery of the decidual source NO (0.27).
  - S-ABS design 2 (v036): real KOLLER_EUR genes.raw kept (z and gene-gene correlations); for each
    replicate the annotation layer is rotated along the genome by an offset uniform in
    [500, n-500] genes (seed 20261107): every marker set's genes and the planted HECA profile move
    together (gene at position i -> i + offset, genome order of genes.raw). Null = rotated sets,
    no planted signal (the real association signal is random with respect to the rotated sets).
    Scenarios: null; dStromal_mid x1, x2, x4; Venous x2; ePV_2 x2; 60 replicates each. The three
    most frequent top sets per scenario and family are printed.
    Readings identical to v0.2.24 (primary 2x; calibrated if null any q < 0.05 <= 0.10 per family;
    ABS relabel; UP recovery >= 0.5), now also for DOWN calibration. If design 2 is calibrated
    and design 1 was not, the v035 calibration failure is attributed to the shifted-z design.
- 2026-11-08 v0.2.26 (after the v036 log; before v037):
  - v036 recorded (calibrated; ABS relabel NO; UP recovery NO at 0.37).
  - H-PAR (v037, rotation design as v036): a lineage-agnostic parenchymal profile PAR = summed HECA
    specificity over stromal and epithelial types (names matching
    Stromal|Fibroblast|SOX9|Lumin|Gland|MUC5B|KRT5|iliated|HOXA13; printed) is planted at 2x and 4x
    (same b_std unit as v036), and a mixture PAR x2 + dStromal_mid x2. 60 replicates each, null rerun.
    H-PAR "supported" if the share of replicates with an endothelial (Arterial / Venous) ABS top set
    under PAR x2 is >= null share + 0.10. Descriptive: share of replicates reproducing the observed
    KOLLER pattern (ABS top endothelial AND UP top stromal: dStromal / eStromal / Fibroblast).
- 2026-11-08 v0.2.27 (after the v037 log): v037 recorded (H-PAR NOT SUPPORTED). The mechanism search
  by single planted sources is closed; no further simulation variants are added post hoc.
  Technical (no rule change): v038 draws draft figures F1-F5 from the v031/v034/v036 outputs
  (rule flip, sign decomposition, joint resolution, simulations, |log2FC| composition).
- 2026-11-08 v0.2.28 (before v039; PLAN v0.4 addendum, docs/plan/PLAN_v0.4_ADDENDUM.md):
  - L1 loci: SNPs with p < 5e-8 in <trait>.pval.txt that are in the MAGMA 1000G EUR reference;
    greedy distance clumping (lead = smallest p; SNPs within 500 kb removed). Candidate genes =
    MAGMA genes (genes.out) whose interval overlaps lead +- 250 kb; prioritized gene = smallest
    MAGMA gene p. Compartments of heca_v014 types by name: stroma (Stromal|Fibroblast), epithelium
    (SOX9|Lumin|Gland|MUC5B|KRT5|iliated|HOXA13), endothelium (Arterial|Venous|Lymphatic),
    perivascular/SM (ePV|mPV|uSMC), immune (uM|uNK|T|B|Mast|DC|Macro|Mono|...); the mapping is
    printed. A gene's compartment share = summed specificity over the compartment's types; the gene
    is assigned to the compartment whose share has the highest percentile across all genes, if
    >= 90th percentile, else "broad". Sensitivity: mean percentile of the top-3 expressed candidates.
    Enrichment: observed locus count per compartment vs the share among all expressed genes
    (one-sided binomial; descriptive, no correction).
    H-MULTI (primary KOLLER_EUR): "supported" if >= 3 compartments each hold >= 10% of the assigned
    loci and no compartment holds > 50%. Reported: whether leiomyoma loci are more concentrated
    (higher top share).
  - L2 leave-one-locus-out (KOLLER_EUR; MAGMA gene-set analysis of the v031 P sets, condition-hide
    ALLGENES): baseline top 3 sets per family (ABS, UP, DOWN) tracked; for each locus all genes
    within lead +- 500 kb are removed from every set except ALLGENES. Driver = locus whose removal
    raises a tracked set's p by >= 10x; then all drivers of a set removed together.
  - L3 (descriptive): a KOLLER_EUR locus is "shared" with leiomyoma / T2D / hypertension if any
    reference SNP within lead +- 250 kb has p < 1e-5 in that trait.
- 2026-11-08 v0.2.29 (after the v039 log; before v040):
  - Technical fix: genes.out chromosome / position columns coerced to numbers (chrX rows made the
    column text, so no gene matched for KOLLER_EUR and the FinnGen traits). L1 rules unchanged.
  - L2b (new): all genome-wide loci of KOLLER_EUR removed at once from every set except ALLGENES;
    reading "signal extends beyond the genome-wide loci" if >= half of the tracked sets keep p < 0.05.
  - L3 baseline (new): 400 random windows of the same width (seed 20261108; centre >= 1 Mb from any
    KOLLER_EUR lead and from chromosome ends; >= 50 reference SNPs); sharing compared by one-sided
    Fisher test per trait.
  - K1b (descriptive): E-MTAB-14039 sample table downloaded and compared with the HECA Mareckova
    libraries (counts by disease / hormone / library type), to reconstruct Koller's 22 / 13 samples.
  - External review (2026-11-08) recorded in docs/plan/PLAN_v0.4_ADDENDUM.md; planned v041 matched
    directional benchmark and v042 external validation will be pre-registered before they run.
- 2026-11-08 v0.2.30 (before v041; matched directional benchmark, answers the external review):
  - Cells and DE: v031 P cells (Mareckova non-hormone, <= 400 per type, seed 20261106), Wilcoxon per
    status x compartment (identity with the v031 thresholded sets checked and printed). DE tables
    and mean CPM per type/status saved (data/processed/genesets/heca_koller_v041_de).
  - Families on the same cells: ABS, UP, DOWN (thresholds as Koller: |log2FC| > 1.2, adj. p < 0.01);
    size-matched ABSN / UPN / DOWNN (N = 400 per type/status, no threshold: top |Wilcoxon score|,
    top score, bottom score); SPEC = continuous specificity (mean CPM share across type/status
    groups) as MAGMA gene-property covariates conditioned on log2 mean expression (direction greater).
  - Design: rotation (as v036): real KOLLER_EUR genes.raw; the annotation layer (sets and covariates)
    and the planted profile rotated together by an offset uniform in [500, n-500] genes, seed 20261108.
    Truths (HECA heca_v014, external to the P cells), strength 2x the R4 unit, targets dStromal_mid,
    Venous, ePV_2: POS = standardized specificity; NEG = its negative (risk in genes depleted in the
    target); EXPR = standardized log2 mean expression in the target (marker-agnostic). Null: 300
    rotations, no planted signal. Truths: 40 rotations each. Replicates stored and resumed.
  - Metrics per family: FWER under the null (any BH q < 0.05 within family) with 95% Clopper-Pearson
    interval, "consistent with <= 5%" if the lower bound <= 0.05; per truth: top-1 exact type (either
    status), top-3 exact, top-1 correct group (stromal / endothelial / perivascular-SM),
    related-state miss (correct group, wrong type), true type q < 0.05.
  - Predictions (stated now): P1 size-matched up >= size-matched |log2FC| (top-1 exact) under POS for
    every target; P2 under NEG, down or |log2FC| (size-matched) beats up for every target; P3 the up
    advantage survives size matching (UPN - ABSN >= 0.2 for >= 2 of 3 targets under POS); P4 no family
    exceeds 0.5 top-1 exact under EXPR for any target (resolution limit).
- 2026-11-08 v0.2.31 (second external review; before v040 and v041 are run):
  - Framing: v040-v041 are prospectively specified follow-up analyses designed after the exploratory
    results v001-v039; pre-specification does not make the programme an independent pre-registration.
  - v040 checks added: autosomes 1-22 only, with non-autosomal SNPs / genes excluded explicitly and
    counted; per GWAS the number of loci with candidate genes and of candidate genes is printed;
    loci are distance-clumped (500 kb, no LD clumping) on GRCh37 positions from the MAGMA 1000G EUR
    reference. L2b wording: "signal persists without the genome-wide loci (distributed over
    sub-threshold regions)" -- not called causal or strictly polygenic. L3 baseline replaced by
    random windows matched on number of genes and reference SNPs (pool of 800 non-overlapping
    windows, 7 nearest controls per locus, seed 20261108); balance printed; Fisher test vs matched
    controls (primary) and vs the whole pool (secondary).
  - v041: P4 (expression-level truths, 0.5 threshold) reclassified as exploratory before the run;
    FWER reported with standard error and 95% interval; differences inside the interval are not
    interpreted. v041 is a conditional benchmark (truths built from expression features); a
    genotype-level simulation independent of expression is out of scope on this hardware and is
    stated as a limitation (Li et al. 2026 used genotype-level simulation).
- 2026-11-08 v0.2.32 (after the v040 log; v042 replaces the unrun second v040 build):
  - v040 recorded. v042 = v040 with the v0.2.31 checks and matched L3 baseline, plus:
  - L2c (formulated after seeing the v040 locus assignments and L2b): loci removed by assigned
    compartment (all loci of a compartment at once). Primary comparisons: endothelium loci vs the
    baseline |log2FC| top set (Venous/CTL); stroma loci vs the baseline up-only top set
    (dStromal_mid/ENDO). Null: 50 random subsets of the same number of loci drawn from the other
    loci (seed 20261108); permutation p = (1 + #null >= observed) / (1 + n). "Supported" if p < 0.05.
  - K1c (descriptive): individuals, (library, individual) pairs and libraries by disease and hormone
    stage, for all libraries, those in the HECA cells file, and the rest.
- 2026-11-08 v0.2.33 (after the v041 log; before v043):
  - v041 recorded (calibrated; P1 YES, P2 YES, P3 NO, P4 exploratory YES).
  - E1: the seven families (v041 definitions, same cells) applied to the real GWAS KOLLER_EUR,
    GCST90205183, ENDO_ALL, LEIOMYOMA, HYPERTENSION, T2D (MAGMA set analysis with condition-hide
    ALLGENES; SPEC as gene-property conditioned on log2 mean expression). Top 3 per family reported.
  - E2 positive controls with expectations fixed now: LEIOMYOMA -> a smooth-muscle / perivascular
    type (uSMCs, mPV, ePV_*) ranked first; HYPERTENSION -> a vascular / perivascular / SM type
    (Arterial, Venous, ePV_*, mPV, uSMCs) ranked first. Reported per family; failure in this
    endometrial atlas is not read as a method error for hypertension (relevant tissues absent).
  - E3 (exploratory): for each v041 scenario, the share of replicates whose top-set groups match the
    observed KOLLER_EUR groups, per family and for all families jointly.
  - Runner: v043 runs the locus atlas v042 first if it has not completed.
- 2026-11-08 v0.2.34 (after the v042 / v043 logs; before v044):
  - v042 / v043 recorded. Koller's cohort = (library, individual) pairs of the HECA cells file with
    hormone-treated samples included (22 / 13; v031 variant H).
  - A (H-STRVEN; formulated after v042 showed the stromal loci weakening the |log2FC| venous set):
    KOLLER_EUR, v031 P sets, MAGMA. Remove the 9 stromal-assigned loci (+- 500 kb) and, separately, the
    8 endothelial-assigned loci; null = 200 random 9-locus subsets of the non-stromal loci (seed
    20261108). "Supported" if for ABS Venous/CTL the stromal removal exceeds the null (permutation
    p < 0.05) AND exceeds the endothelial removal. Also reported for ABS Arterial/CTL, Venous/ENDO,
    DOWN Ciliated/ENDO, Luminal/ENDO and UP dStromal_mid/ENDO.
  - B external atlas outside the reproductive tract: GTEx v9 snRNA-seq (Eraslan et al. 2022; CELLxGENE
    collection a3ffde6c-7ad2-498a-903c-d58e732f7470; largest dataset). Cell types = CELLxGENE
    cell_type with >= 50 cells and >= 3 donors; <= 400 cells per type (seed 20261108); Wilcoxon each
    type vs the rest (one compartment); families as v041 (ABS / UP / DOWN thresholds, size-matched
    N = 400, continuous specificity). Six GWAS as v043. Positive controls fixed now: LEIOMYOMA ->
    a type matching "smooth muscle" ranked first; HYPERTENSION -> a type matching
    "smooth muscle|endothel|pericyte|vascular" ranked first. T2D -> adipocyte is exploratory.
- 2026-11-08 v0.2.35 (after the v044 log; before v045):
  - v044 recorded. Wording fix: the GTEx v9 atlas lacks uterine / female reproductive-tract tissue but
    includes breast and prostate; it is described as "an atlas without uterine or female
    reproductive-tract tissue", not "non-reproductive".
  - C1 (post hoc label check): identity of ambiguous GTEx types ("contractile cell", "smooth muscle
    cell", "enteric smooth muscle cell", "fibroblast", "basal-myoepithelial cell of mammary gland",
    "professional antigen presenting cell") from author annotation columns, tissue and canonical
    marker panels. "contractile cell" counts as vascular mural only if the share of mural-matching
    author labels (pericyte|SMC|smooth|vascular|mural) is >= 0.5 in every author annotation column;
    marker rule only if no such column exists. The pre-registered v044 hypertension verdict (SPEC only)
    remains the primary result; the post hoc verdict is reported separately.
  - C2 H-NONSPEC (formulated after v044 showed the same down-set cell types first for every trait):
    traits KOLLER_EUR, LEIOMYOMA, HYPERTENSION, T2D (6 pairs); per family, Spearman correlation of
    -log10 p across the labels common to all families (v043 HECA, v044 GTEx rankings), averaged over
    pairs. Supported if mean(DOWN, DOWNN) - mean(UP, UPN, SPEC) > 0 with a bootstrap 95% CI above 0
    (1,000 label resamples, seed 20261108) in both atlases.
  - C3 H-ABSDOWN: supported if the median depleted-gene share of GTEx |log2FC| sets is > 0.5 AND, in
    both atlases, the Spearman correlation of |log2FC| rankings with DOWN rankings exceeds that with UP
    rankings in >= 4 of 6 traits.
  - C4 H-GENERIC: MAGMA gene-property analysis (direction=greater) with binary set indicators for all
    six set families plus the specificity shares, both atlases, six GWAS. Universe: genes with CPM >= 1
    in >= 1 group and with gnomAD v2.1.1 LOEUF (if the download fails: expression covariates only).
    Base = no conditioning (sets) / Average (SPEC); conditioned = Average (log2 mean CPM) + Breadth
    (share of groups with CPM >= 1) + LOEUF. Supported if the pooled share of BH-significant labels
    retained after conditioning is lower for down-type (DOWN, DOWNN) than up-type (UP, UPN, SPEC)
    families in both atlases. H-DSTROMA: KOLLER_EUR UP and UPN dStromal_mid/ENDO keep p < 0.05 after
    conditioning (HECA). Also reported: positive controls and key sets (ABS Venous/CTL, Arterial/CTL;
    DOWN Luminal/ENDO, Ciliated/ENDO) before / after.
- 2026-11-08 v0.2.36 (after the v045 log; before v046):
  - v045 recorded (H-NONSPEC, H-ABSDOWN, H-GENERIC not supported; H-DSTROMA supported).
  - Deviation (C1 rule): the fixed rule required a mural majority in every author annotation column,
    including the coarse hierarchy columns "Cell types level 2/3" (Muscle / Stromal), which never use
    mural vocabulary; it therefore printed "not vascular mural". The specific annotations (Broad cell
    type Pericyte/SMC 98%, Granular Pericyte/SMC I/II 100%) and the marker panels identify "contractile
    cell" as Pericyte/SMC. From now on it is reported as Pericyte/SMC. The pre-registered v044
    hypertension verdict (SPEC only) stays primary; the corrected-label verdict (vascular first with
    ABS, UP, ABSN, UPN, SPEC; not DOWN, DOWNN) is post hoc.
  - Correction: "basal-myoepithelial cell of mammary gland" is 71% oesophageal-mucosa myofibroblasts;
    the v044 reading "hormone-responsive mammary epithelium" is withdrawn.
  - D (formulated after v045; K7 rule from v0.2.23, MAGMA instead of S-LDSC): universe as v045 C4
    (CPM >= 1 in >= 1 group, gnomAD LOEUF available); |log2FC| sets split into the UP part and the DOWN
    part (= the UP and DOWN family sets); MAGMA gene-property analysis of binary indicators,
    direction=greater.
    - D1: per trait, Spearman of the depleted-gene share with -log10 p of the |log2FC| set (descriptive).
    - D2 H-K7REP: for the top-5 |log2FC| sets per trait (and the KOLLER_EUR key sets in HECA), joint
      analysis UP part | DOWN part and DOWN part | UP part; class by K7 (p < 0.05). Supported if, among
      top-5 sets with |log2FC| p < 0.05 pooled over the 6 traits, >= 50% are DOWN-carried (or DOWN-only)
      in both atlases.
    - D3 H-MASK: leiomyoma (pre-registered smooth-muscle patterns): supported if, in both atlases, the
      expected group is not first under |log2FC| and >= half of the |log2FC| sets ranked above it (at
      most 10) are DOWN-carried. Hypertension is reported descriptively (GTEx with the corrected
      Pericyte/SMC label).
- 2026-11-08 v0.2.37 (after the v046 log; before v047):
  - v046 recorded (H-K7REP and H-MASK not supported; masking of the true source replicated in MAGMA).
  - X (exact Koller cohort): Koller's 22 / 13 samples = v031 variant H (K1c, v042). So far the rule
    comparison, sign decomposition and joint resolution in S-LDSC used variant P; they are repeated on
    variant H with the v031 H marker sets (|log2FC| LD scores reused; up-only and down-only LD scores
    new). Traits KOLLER_EUR, GCST90205183, ENDO_ALL, LEIOMYOMA, HYPERTENSION (T2D omitted for run time).
  - H-EXACT supported if (a) the first KOLLER_EUR |log2FC| annotation is endothelial (Venous or
    Arterial), (b) the first KOLLER_EUR up-only annotation is decidual stroma (dStromal_*), and (c) for
    leiomyoma the smooth-muscle group (uSMCs or mPV, pattern of v0.2.21) is first with up-only and not
    first with |log2FC|.
  - K2 / H-ABS, K4 (KOLLER_EUR, GCST90205183), K6 (rules v0.2.21) and K7 H-SIGN, K8 (rules v0.2.23;
    traits KOLLER_EUR, GCST90205183) are repeated unchanged on variant H.
- 2026-11-08 v0.2.38 (after the v047 log; before v048):
  - v047 recorded (H-EXACT supported; K4 / K7 / K8 on variant H as on variant P; H-ABS and H-SIGN
    all-of rules not met on variant H).
  - v048 draws manuscript figures from finished runs only (no new analyses): F1 rule flip (exact
    cohort, S-LDSC), F2 sign decomposition / composition / masking, F3 joint resolution, F4 benchmark,
    F5 positive controls (S-LDSC exact cohort; MAGMA HECA; MAGMA GTEx v9 with the pre-registered and
    the post hoc Pericyte/SMC labels shown as separate rows), F6 stromal loci, S1 conditioning.
    Colours: categorical slots 1-3 (validated all-pairs, light surface) for compartment groups;
    sequential blue for heat maps.
