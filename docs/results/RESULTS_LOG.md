# Paper 9 — Results log

Each entry: date, run, hypothesis, pre-registered rule, result, verdict. Entries are never
edited after the fact; corrections are added as new entries.

## 2026-11-05 · run p9_v004_ldsc_h2_rg_20261105T210005Z · G3 part A (rg = 1 between subtypes)

Tool: LDSC Python 3 port (belowlab 2-to-3), upstream nose suite 167/167 PASS.
Data: FinnGen R13, HapMap3 no-MHC (1,190,321 SNPs after munge).

SNP-heritability (observed scale; all subtypes pass the h2 z >= 4 gate):

| Trait | Cases | h2_obs (SE) | z | Intercept (SE) |
|---|---|---|---|---|
| ENDO_OVARY | 8,002 | 0.0524 (0.0060) | 8.7 | 1.113 (0.011) |
| ENDO_PERITONEUM | 7,758 | 0.0428 (0.0053) | 8.1 | 1.083 (0.009) |
| ENDO_DEEP | 3,890 | 0.0123 (0.0017) | 7.2 | 1.037 (0.009) |
| ENDO_RECTOVAG | 3,263 | 0.0283 (0.0040) | 7.1 | 1.038 (0.009) |
| ENDO_STAGE34 | 10,284 | 0.0243 (0.0028) | 8.7 | 1.106 (0.012) |
| ADENOMYOSIS | 5,864 | 0.0250 (0.0043) | 5.8 | 1.056 (0.009) |
| LEIOMYOMA | 42,816 | 0.0487 (0.0058) | 8.4 | 1.164 (0.016) |
| HYPERTENSION | 162,639 | 0.0740 (0.0043) | 17.2 | 1.496 (0.025) |
| T2D | 89,727 | 0.0749 (0.0044) | 17.0 | 1.482 (0.026) |
| GCST90205183 | (N=474,160) | 0.0123 (0.0012) | 10.3 | 0.953 (0.007) |

ENDO_ALL liability-scale h2 = 0.150 (K = 0.10); observed-scale value not parsed (LDSC prints
liability only when prevalences are given) — to be re-read from the LDSC log.

G3 primary family (two-sided Wald test of rg = 1, Holm over 6 pairs):

| Pair | rg [95% CI] | p | Holm |
|---|---|---|---|
| OVARY ~ PERITONEUM | 0.893 [0.807, 0.978] | 0.014 | 0.084 |
| OVARY ~ DEEP | 0.880 [0.741, 1.020] | 0.093 | 0.454 |
| OVARY ~ RECTOVAG | 0.884 [0.749, 1.018] | 0.091 | 0.454 |
| PERITONEUM ~ DEEP | 1.011 [0.861, 1.160] | 0.886 | 1 |
| PERITONEUM ~ RECTOVAG | 1.032 [0.880, 1.183] | 0.682 | 1 |
| DEEP ~ RECTOVAG | 0.994 [0.973, 1.014] | 0.552 | 1 |

Verdict (pre-registered rule): **G3 part A NOT SUPPORTED** — rg = 1 not rejected for any pair.
Pattern (descriptive): the ovarian subtype is the only one with rg < 1 point estimates
(0.88–0.89 vs every other subtype); peritoneal, deep and rectovaginal are ~1.
Caveats: (i) DEEP ~ RECTOVAG is nested (rectovaginal is a subset of deep; gcov_int 0.94) and
uninformative — dropping it leaves the minimum Holm p at 0.07; (ii) case overlap across subtype
codes biases rg towards 1, so the test is conservative for "distinct disease" claims.

Descriptive rg (G5 context): ENDO_ALL ~ ADENOMYOSIS 0.816 (0.045); ~ LEIOMYOMA 0.495 (0.043);
~ HYPERTENSION 0.163 (0.040); ~ T2D 0.090 (0.039); FinnGen ENDO_ALL ~ GCST90205183 0.890 (0.050),
gcov_int 0.051. Leiomyoma correlates more with ovarian/peritoneal (0.41/0.46) than deep/rectovaginal
(0.26/0.29). Deep and rectovaginal correlate less with GCST90205183 (0.76/0.73) than ovarian/peritoneal
(0.88/0.92).

Next for G3: part B (subtype-specific cell-context enrichment) in the G1 pipeline; exploratory X2.

## 2026-11-06 · run p9_v012_ldsc_seg_20261106T154700Z · G4 check 1 (LDSC-SEG, GSE266265 CTL sets)

Tool checks: partitioned LD scores for 59 sets (22 chromosomes); M_5_50 rebuilt and matched the
baseline v1.2 count exactly on chr22 (83,323); MAGMA sanity PASS (WNT4, GREB1, KDR, SYNE1 ranks
1-4; 11 known endometriosis genes Bonferroni-significant).

G4 primary (GCST90205183, A_refined_CTL, one-sided, Holm over 4):

| Set | coef (SE) | p | Holm |
|---|---|---|---|
| dS4 (5 CTL donors) | 1.71e-9 (7.90e-10) | 0.015 | 0.061 |
| dS3 (12) | 1.28e-9 (6.71e-10) | 0.029 | 0.085 |
| dS2 (21) | -1.91e-10 (6.01e-10) | 0.625 | 1 |
| dS1-myofibroblast (22) | -1.02e-9 (5.76e-10) | 0.962 | 1 |

Verdict (pre-registered rule): **check 1 NOT REPLICATED** at Holm < 0.05.
Interpretation limits: dS4 and dS3 are ranks 1-2 of 42 types (direction consistent with HECA);
no type in any trait passes BH (best GCST q = 0.60); the hypertension positive control ranks
vascular smooth-muscle / perivascular types 1-4 (Prv-STEAP4, VSMC, fib C7 secretory, Prv_VSMC
secretory) but q = 0.26 — the within-tissue LDSC-SEG design is power-limited. Post hoc (not
pre-registered): P(two decidual types occupy ranks 1-2 of 42 by chance) = 6/861 = 0.007.

Descriptive tops (no set passes BH): ENDO_OVARY dS4 (p = 0.0034); ENDO_PERITONEUM LEC, EC-artery,
EC-tip; ENDO_DEEP mesothelial, LEC; LEIOMYOMA and T2D no signal (myometrium / pancreas absent).

## 2026-11-06 · run p9_v013_magma_property_20261106T162822Z · G4 check 2 (MAGMA gene-property)

Continuous specificity (A_refined_CTL, 14,768 genes, 42 types), conditioned on average expression,
one-sided. GCST90205183 decidual types (Holm over 4):

| Set | beta (SE) | p | Holm | rank / 42 |
|---|---|---|---|---|
| dS3 | 0.824 (0.501) | 0.050 | 0.201 | 3 |
| dS1-myofibroblast | 0.683 (0.560) | 0.111 | 0.334 | 5 |
| dS4 | 0.382 (0.401) | 0.171 | 0.341 | 6 |
| dS2 | 0.185 (0.527) | 0.362 | 0.362 | 19 |

Verdict: **check 2 NOT REPLICATED.** Best GCST types overall: T_N/T_CM (p 0.013), T_EM (0.028);
nothing passes BH (best q 0.53).
Positive control PASSED: hypertension (FinnGen, N = 500,103) Prv_VSMC secretory_late q = 0.006,
EC-artery 0.029, Prv_VSMC secretory 0.029, EC-tip 0.033, VSMC 0.033 (4 of top 5 vascular);
family B stromal@secretory q = 0.004, endothelial@secretory q = 0.007.
Interpretation: the method detects concentrated cell-type heritability when it exists; endometriosis
heritability shows no comparable concentration in any eutopic cell type of GSE266265. Decidual types
remain in the upper part of the ranking in both methods (consistent direction, weak signal).

FinnGen, exploratory (no set passes BH): ENDO_OVARY VSMC, fib C7 secretory_late, EC-artery, dS4
(q 0.27); ENDO_PERITONEUM NK3, LEC, NK1, cDC1; ENDO_DEEP dS2 (p 0.0049, q 0.20), eF3, eF4-CXCL14,
eF2; ENDO_ALL none; LEIOMYOMA none (no myometrium in GSE266265). H-sub is not consistent across
LDSC-SEG and MAGMA.

## 2026-11-06 · run p9_v014_heca_magma_20261106T192726Z · G4 check 3 (HECA cell types, MAGMA)

HECA control donors without hormones: 30 donors, 159,904 cells (136,670 superficial + 23,234
whole uterus), 32 eligible types, universe 14,348 genes.

GCST90205183, Holm over HECA decidual states:

| Set | beta (SE) | p | Holm | rank / 32 |
|---|---|---|---|---|
| dStromal_mid | 0.580 (0.416) | 0.082 | 0.245 | 1 |
| dStromal_late | 0.491 (0.405) | 0.113 | 0.245 | 3 |
| dStromal_early | 0.446 (0.369) | 0.113 | 0.245 | 4 |

Verdict: **check 3 NOT REPRODUCED.**
Positive controls on the same atlas and method: LEIOMYOMA -> uSMCs rank 1, beta 0.944 (0.238),
p 3.7e-5, q 0.001 (FinnGen h2 z = 8.4, i.e. weaker than GCST h2 z = 10.3); HYPERTENSION -> mPV
(q < 0.001), Arterial (q < 0.001), ePV_1b (q < 0.001), Lymphatic (q 0.021), ePV_2 (q 0.034).

Synthesis of G4 checks 1-3: decidual stromal programs rank at or near the top for endometriosis in
two atlases and two methods but never reach significance, whereas the identical pipeline detects
strong cell-type concentration for leiomyoma and hypertension. Endometriosis SNP-heritability is not
concentrated in any eutopic endometrial cell type at the level these methods detect.

Exploratory: ENDO_OVARY -> ePV_1b (p 0.0032, q 0.10), uSMCs (q 0.12), Arterial (q 0.12), consistent
with GSE266265 MAGMA (VSMC, EC-artery, Prv_VSMC) -> X3, with leiomyoma comorbidity as a candidate
confounder (rg OVARY-LEIOMYOMA 0.41).

## 2026-11-06 · run p9_v016_gse213216_prepare_20261106T193904Z · G4 check 4 + G1 data QC

HECA fGWAS input (from the authors' repository): GCST90205183 (GRCh38), per-cell-type mean
expression at gene TSS -> same GWAS as ours.

Check 4a: the 27 Bonferroni MAGMA genes in the HECA universe are not dStromal-specific
(dStromal_late rank 4/32, p_emp 0.20; early 14/32; mid 23/32; best type SOX9_functionalis_I p_emp
0.035, one of 32).
Check 4b: removing 280 genes within +-1 Mb of those loci: dStromal_mid beta 0.580 -> 0.501
(rank 1 -> 2), early 0.446 -> 0.347 (4 -> 6), late 0.491 -> 0.394 (3 -> 5).
Interpretation: the weak decidual tilt is polygenic, not driven by a few strong loci; this
explanation of the HECA (fGWAS) vs polygenic-method discrepancy is not supported.

GSE213216 preparation: 46 scRNA samples processed (2 Visium excluded); 828 clusters; lineage
agreement 90.6%; fine-label margin median 0.012. Patients with >= 10 cells per lineage: eutopic 9
(all lineages); peritoneal lesion 12 (epithelial 7); endometrioma 5-6 (epithelial 1).

## 2026-11-06 · run p9_v017_g1_20261106T194734Z · G4 sensitivity (HECA all donors)

HECA all 63 donors (any pathology, hormones included), 33 eligible types, MAGMA gene-property on
GCST90205183: dStromal_mid rank 1/33 (beta 0.640, p 0.081), dStromal_early rank 5 (p 0.16),
dStromal_late rank 6 (p 0.20). Donor selection does not explain the HECA (fGWAS) vs polygenic
discrepancy. (G1 statistics of this run were not produced: "less" direction unsupported by MAGMA.)

## 2026-11-06 · run p9_v018_g1_20261106T195443Z · G1 primary + secondary

Rule (CHANGELOG v0.2.7 + v0.2.9, v018 technical): peritoneal lesion vs eutopic endometrium,
5 lineages, GCST90205183, MAGMA gene-property (D = lesion - eutopic log2CPM, conditioned on AVG),
Holm over lineages; supported if Holm p_greater < 0.05 with beta > 0; alternative if Holm
p_less < 0.05. Patients lesion/eutopic: stromal 12/9, epithelial 5/9, endothelial 10/6,
myeloid 10/8, lymphocyte 10/5; 11,246-12,630 genes.

Primary: no lineage shows lesion-up enrichment (p_greater 0.50-0.997, Holm 1.0).
Eutopic-up: lymphocyte beta -0.055 (se 0.020) p_less 0.0035 Holm 0.017; epithelial Holm 0.063;
stromal 0.101; myeloid 0.101; endothelial 0.504.
Verdict by rule: G1 ALTERNATIVE (lymphocyte). The lesion-niche hypothesis is not supported.

Secondary (descriptive), lesion vs eutopic:
- HYPERTENSION: same direction, stromal Holm 0.028, epithelial 0.028, myeloid 0.011, lymphocyte
  0.028 (beta -0.062 vs -0.055 for endometriosis); endothelial 0.31.
- LEIOMYOMA: no lineage below Holm 0.22; lymphocyte beta +0.014.
- ENDO_PERITONEUM (subtype matched to the lesions): nothing (min Holm 0.055, endothelial, negative).
- ENDO_ALL: epithelial Holm 0.006 (negative); ENDO_OVARY: epithelial 0.005, myeloid 0.038 (negative).
Endometrioma vs eutopic: GCST90205183 nothing (min Holm 0.32); ENDO_OVARY nothing; HYPERTENSION
stromal 0.017, myeloid 0.021, lymphocyte 0.021 (negative). Paired patients (4, stromal only):
beta -0.022, p_less 0.052.

Interpretation: the only signal is the pre-specified alternative direction, and it is shared by an
unrelated trait (hypertension) at equal or larger effect, absent for the subtype-matched GWAS and
for the endometrioma contrast; it is therefore not interpreted as endometriosis-specific (post hoc;
diagnostics G1-S in v019). No evidence that endometriosis heritability is anchored in
lesion-context programs.

## 2026-11-06 · run p9_v019_tissue_g1s_20261106T201024Z · G6 (exploratory) + G1-S (post hoc)

G6, GTEx v8, 52 tissues, MAGMA gene-property | Average, 16,796 genes (rule v0.2.10; re-classified
exploratory in v0.2.11 before this log was seen).
GCST90205183: Uterus p 3.7e-07, Fallopian Tube 3.6e-06, Cervix-Ectocervix 1.9e-05,
Cervix-Endocervix 2.4e-05, Ovary 3.2e-03 (q 0.033), Vagina 0.025; top 5 all reproductive.
Verdict by rule: REPRODUCTIVE (all 6 reproductive tissues Holm < 0.05).
ENDO_ALL same pattern; ENDO_OVARY strongest (Fallopian Tube 1.1e-09, Uterus 2.4e-09);
ENDO_PERITONEUM Fallopian Tube / Uterus only; ENDO_DEEP weaker (Fallopian Tube 1.9e-03,
Nerve-Tibial 2.5e-03, Colon-Sigmoid 6.1e-03; q 0.065-0.089).
Controls: LEIOMYOMA Uterus rank 1 (PASS) with Esophagus-Muscularis rank 3 (smooth muscle);
HYPERTENSION top = Uterus (p 0.0056, q 0.135), heart tissues; artery not in top 3
(pre-specified positive control FAILED); T2D: Brain-Cerebellum (q 0.025), pancreas rank 4;
reproductive tissues not enriched.
Reading: at tissue level the polygenic endometriosis signal is concentrated in reproductive-tract
tissues, broader than for leiomyoma or hypertension (cervix, ovary, vagina significant only for
endometriosis); calibration caveat: the hypertension control failed.

G1-S (post hoc, descriptive):
S1 T2D (lesion vs eutopic): epithelial beta -0.049 Holm 0.001; stromal p_less 0.040; lymphocyte
beta +0.010 (no tilt). The eutopic-up tilt occurs for T2D in epithelium, for HYPERTENSION in
4/5 lineages, for GCST90205183 in lymphocytes (others nominal), not for LEIOMYOMA.
S2 LOEUF + tau conditioning: lymphocyte GCST beta -0.048 -> -0.045 (6% shrink), HYPERTENSION
-0.063 -> -0.055 (13%) -> generic gene-property explanation NOT SUPPORTED (fixed criterion).
S3 GCST lymphocyte conditional on HYPERTENSION gene z: beta -0.055 -> -0.051 (Holm 0.031);
on T2D gene z: -0.053 (Holm 0.023).
S4 lymphocyte composition: eutopic NK2+NK3 23%, T_N/T_CM 35%; peritoneal lesion Treg 18%,
NK 9%; endometrioma T_EM 74%. D for lymphocytes therefore partly encodes NK / naive-T vs Treg /
effector-memory composition.
Reading: the lymphocyte tilt is not explained by constraint, expression breadth or shared gene
z with hypertension / T2D, but it is shared in direction and size with hypertension (a known
comorbidity, rg 0.16) and similar eutopic-up tilts occur for unrelated traits in other
lineages; it remains "not shown to be endometriosis-specific" (A7 difference test pending).

## 2026-11-06 · run p9_v020_g4c5_limits_20261106T205328Z · Fonseca check, G4 check 5, detection limits

Fonseca table (authors' GEO_id_patient_class.csv): 46/46 of our GSMs present; patient grouping
identical (0 conflicts); classes agree with our contexts. 3 authors' GSMs (GSM6634272-4: patient 2
endometrioma x2, patient 3 eutopic) were never processed by us (absent from the GEO RAW tar).
No per-cell annotation files in the repository (Seurat objects on Box only).

G4 check 5 (HECA control donors, 32 types). GCST90205183: M1 specificity|Average dStromal_mid
rank 1 p 0.082 (Holm 0.25); M2 log-mean|Average top eStromal_MMPs p 0.068; M3 log-mean without
conditioning: all three dStromal types Holm 0.006 (eStromal_MMPs, dStromal_late, uSMCs top-3).
Under M3 the dStromal types are also Holm-significant for HYPERTENSION (0.017), LEIOMYOMA (0.0005)
and T2D (0.0008). Verdict by the fixed rule: NOT SUPPORTED (the rule required >= 2 shared top-3
types with hypertension; shared = 0, hypertension M3 top-3 are vascular). Post hoc observation:
an unconditioned mean-expression annotation makes decidual stroma "significant" for every trait
tested, i.e. significance under M3 is not trait-specific. Also M1: T2D dStromal_late Holm 0.019.

Detection limits (80% power, one-sided alpha 0.05/k): G4 dStromal MDE BETA_STD 0.025; observed
dStromal_mid BETA_STD 0.012 (upper 95% bound ~0.029); positive controls in their own GWAS:
LEIOMYOMA -> uSMCs BETA_STD 0.035, HYPERTENSION -> mPV 0.045. A decidual effect of the size of
the positive controls is therefore outside the 95% CI. G1 MDE BETA_STD 0.027-0.030 per lineage;
lymphocyte observed -0.026.

## 2026-11-06 · run p9_v021_audit_20261106T205327Z · audit A0-A7

A0 Koller 2026 Zenodo (open, CC-BY-NC-4.0): endometriosis_EUR_wo_23andMe (169.7 MB),
endometriosis_wo.adenomyosis_EUR, adenomyosis_EUR, endometriosis_EAS / AFR / AMR, PRS weights.
A1 munge writes all 1,190,321 merge-allele SNPs; non-missing 1,153,655-1,175,329 per trait
(GCST 1,175,329 of 1,183,085 pre-filtered) -> the v004 "SNPs" number was the merge-list size
(logging artifact, no data error).
A2 alleles: same A1 in 100% of shared SNPs; GWS GCST SNPs: sign concordance 1.000 with all FinnGen
endometriosis endpoints (n 137-190), 0.90 with leiomyoma, 0.55 hypertension, 0.25 T2D (n 32);
top GCST SNPs identical p in munged and MAGMA inputs. MAGMA used rsID synonyms (6.0 M).
A3 LDSC: FinnGen ratio 0.32-0.41 for all endpoints (likely Finnish-LD vs 1000G EUR mismatch and
relatedness; rg uses free intercepts); GCST90205183 intercept 0.953, mean chi2 1.09 (consistent
with genomic-control-corrected meta-analysis statistics: conservative, lower power).
A4 GSE266265 counts: 92.4% integer, fractional mass concentrated at .5, min 1e-08; column sums
equal nCount_RNA exactly -> deposited counts slot is ambient-corrected-type non-integer data from
which the authors' totals were computed; CPM-based specificity unaffected (no count-model DE used).
A5 GSE266265 leave-one-control-donor-out (25 donors): no dS / stromal type reaches Holm < 0.05 in
any run; beta > 0 in 25/25 for all; dS3 p 0.030-0.084.
A6 G1 lymphocyte: no lineage filter Holm 0.03; >= 50 cells Holm 0.02; leave-one-unit-out p_less
< 0.05 in 15/15; patient-label permutation (B = 100) p_less 0.020 vs parametric 0.0035 (permutation
t SD 1.37: parametric MAGMA p-values for donor-estimated contrasts are anti-conservative);
HYPERTENSION permutation p 0.050.
A7 cross-trait difference (rho = LDSC gcov_int 0.001-0.003): GCST vs HYPERTENSION no lineage
differs (|z| <= 1.07; lymphocyte z 0.22); GCST vs T2D lymphocyte z -2.01, p 0.045 (uncorrected).

## 2026-11-06 · run p9_v022_x4_20261106T205327Z · X4 organ identity (exploratory)

X4a MULL | Average + SMOOTH: GCST90205183 beta 0.067 p 1.2e-05 (Uterus alone 3.3e-06);
ENDO_OVARY 4.6e-06; ENDO_ALL 3.1e-04; ENDO_PERITONEUM 0.013; ENDO_DEEP 0.029; LEIOMYOMA 0.006;
HYPERTENSION 0.020; T2D 0.21 (but T2D Uterus | SMOOTH 0.0014).
X4b by HECA breadth tertile (GCST): T1 broad p 0.094, T2 3.9e-04, T3 specific 0.056, T1 > T3
p 0.34; FinnGen endometriosis traits strongest in T3.
X4c GCST dStromal_mid p 0.079 -> 0.051 after conditioning on MULL (beta rises).
Verdict by the fixed rule: Muellerian beyond smooth muscle YES; endometriosis-specific NO
(hypertension p 0.020); organ-identity pattern NO (falsified as formulated).

## 2026-11-06 · run p9_v023_koller_20261106T211301Z · larger, partially overlapping GWAS (Koller 2026 EUR without 23andMe)

Data: Zenodo endometriosis_EUR_wo_23andMe (UKB, AoU, MVP, FinnGen, EstBB, IEGC; METAL sample-size
weighted; columns SNP, Allele1, Allele2, Freq1, N, Z, P); 8,779,723 rsID SNPs; per-SNP N ~118,663
(sum of study weights). LDSC: h2_obs 0.090, lambda GC 0.972, mean chi2 1.095, intercept 0.80
(GC-corrected, conservative). rg with GCST90205183 0.958 (se 0.021; gcov_int 0.53 = heavy sample
overlap), with FinnGen ENDO_ALL 0.957 (gcov_int 0.69). MAGMA: 72 Bonferroni genes (SYNE1, GREB1,
RNLS, WNT4, CDC42, ESR1, KDR, RSPO3, FSHB ...); 27/28 GCST Bonferroni genes also Bonferroni here;
gene Z correlation 0.505.

Re-runs (rules v0.2.14: holds = same sign and p < 0.05, Holm where the original used Holm):
- G6: Uterus p 1.1e-12, Fallopian Tube 2.5e-10, Cervix-Ecto 1.4e-08, Cervix-Endo 3.8e-08, Ovary
  4.5e-05, Vagina Holm 6.3e-04 -> HOLDS (stronger).
- X4a MULL | Average + SMOOTH: beta 0.101, p 9.3e-10 (Uterus 3.9e-11) -> HOLDS (stronger).
  X4b: T1 p 0.012, T2 1.1e-05, T3 3.2e-05; BETA_STD 0.186 / 0.222 / 0.179 -> still no breadth
  gradient (organ-identity pattern still NO).
- G4 HECA M1 (specificity | Average): dStromal_mid rank 1, p 0.012, Holm (3 dStromal) 0.037 ->
  by the fixed rule G4 is REPRODUCED with the larger GWAS (dStromal_early 0.079, late 0.14);
  ranks 2-3 are ePV_2 (p 0.014) and Arterial (0.017); across 32 types not significant (BH).
  M3 (no conditioning): dStromal Holm ~1e-05 (non-specific under M3, as for all traits in v020).
  Caveat from v020: T2D also reaches Holm < 0.05 for dStromal under M1 (dStromal_late 0.019,
  dStromal_mid 0.044).
- G1 primary: lymphocyte beta -0.025, p_less 0.12 (Holm 0.36; permutation p 0.22) -> the GCST
  lymphocyte alternative does NOT hold; epithelial beta -0.044, Holm < 0.01 (as T2D, ENDO_ALL,
  hypertension in earlier runs). No lineage shows lesion-up enrichment (Holm_greater 1.0).
Reading: the tissue-level reproductive-tract signal is robust to the GWAS; the decidual cell-type
signal becomes nominally/Holm-significant only with the larger GWAS and has trait-specificity
caveats; eutopic-up lineage tilts move between lineages across GWAS and occur for unrelated traits.

## 2026-11-06 · run p9_v024_g4s_20261106T214033Z · G4-S decidual specificity + G1 epithelial permutation

G4-S (HECA control donors, M1, 14,148 genes with all covariates), dStromal_mid:
KOLLER_EUR base p 0.016; | + T2D gene z 0.026; | + HYPERTENSION z 0.021; | + LEIOMYOMA z 0.023;
| + ePV_2, Arterial, uSMCs specificity 0.109 (beta 0.95 -> 0.63); | + MULL 0.027.
GCST90205183 base 0.067; + MULL 0.045; others 0.075-0.086.
Readings (fixed in v0.2.15): KOLLER - not explained by generic gene-level signal YES; separable
from perivascular / smooth-muscle programs NO; beyond Muellerian expression YES.
M1 panel (top type; dStromal Holm over 3): KOLLER dStromal_mid (Holm 0.037); GCST dStromal_mid
(0.245); ENDO_ALL ePV_2; ENDO_OVARY ePV_1b; ENDO_PERITONEUM Arterial; ENDO_DEEP sHormones;
LEIOMYOMA uSMCs; HYPERTENSION mPV; T2D ePV_1b with dStromal_late Holm 0.019 and dStromal_mid
Holm 0.044 -> the G4 criterion is also met by T2D.
Reading: within eutopic endometrium the endometriosis signal leans to a mesenchymal-perivascular
continuum (decidual stroma, perivascular, arterial / smooth muscle) rather than to decidual stroma
specifically; the decidual "hit" is not separable from perivascular programs and an unrelated
trait meets the same criterion -> not endometriosis-specific.

G1-S epithelial contrast (5 lesion vs 9 eutopic units), permutation B = 100:
KOLLER t -3.63, permutation p 0.069 (parametric 1.5e-04; permutation t SD 2.08); GCST t -2.15,
p 0.109 (parametric 0.016; SD 1.47); T2D t -3.55, p 0.040 (parametric 1.9e-04; SD 1.73).
Reading: parametric MAGMA p-values for donor-estimated contrasts can be anti-conservative by
two to three orders of magnitude; after calibration no epithelial endometriosis tilt remains.

## 2026-11-06 · runs p9_v026_env_scdrs_20261106T220851Z + p9_v025_scdrs_20261106T221011Z · scDRS (Ochoa 2025 claim)

Environment: main .venv intact (numpy 2.0.2; all pipeline modules import); only scdrs failed there
(scikit-misc cp39 wheel built for numpy 1.x). scDRS ran in .venv_scdrs (numpy 1.26.4, scanpy 1.10.3).
Cells: per-sample re-processing (v016 rules), lineage-agreement filter, <= 150 cells per sample x
label: 42,891 cells (peritoneal lesion 24,847, endometrioma 7,500, eutopic 10,544).
Lesion (peritoneal + endometrioma, 32,347 cells), scDRS association MC p (1,000 control sets):
- Macrophage labels (MPhi4-infiltrated, MPhi5-activated) + myeloid lineage: GCST90205183 min Holm
  0.665; KOLLER_EUR 0.977; HYPERTENSION myeloid p 0.003 (Holm 0.009), MPhi5 0.012; T2D 0.15;
  LEIOMYOMA 1.0. Verdict by the fixed rule: Ochoa claim NOT reproduced; the only myeloid signal is
  for hypertension.
- Lineages: GCST epithelial 0.033, stromal 0.042; KOLLER stromal 0.023, epithelial 0.020,
  endothelial 0.038; LEIOMYOMA stromal 0.038; lymphocytes p > 0.95 for every trait.
- Top labels: GCST eF2 0.001, lumenal 1 0.012; KOLLER dS4 0.001, eF2 0.001, dS2 0.009;
  T2D eF2 0.007; LEIOMYOMA eF2 0.001, dS2 0.006; HYPERTENSION pre-cDC2 0.001, EC-capillary 0.003.
Eutopic (10,544 cells): stromal lineage GCST 0.032, KOLLER 0.037, LEIOMYOMA 0.009; KOLLER fib C7
0.001; HYPERTENSION myeloid 0.038.
Reading: with two larger endometriosis GWAS the lesion-macrophage prioritization does not reproduce;
the endometriosis scores lean to endometrial-type stromal states (eF2, dS2, dS4), which leiomyoma
and T2D share -> again not trait-specific. Caveats: our labels and sample definitions differ from
Ochoa's (they re-classified 4 peritoneum samples and re-annotated immune clusters); cells, not
donors, were the units; no correction across labels -> Ochoa-faithful re-run with donor-level
tests registered (v0.2.18).

## 2026-11-06 · run p9_v027_resolution_20261106T222520Z · Option A core (R1-R5, B0)

R1 (HECA control-donor specificity): dStromal_mid correlates with dStromal_early / late (Pearson
0.84 / 0.68), ePV_2 (0.35; Spearman 0.67), uSMCs (0.20 / 0.54), not with Arterial (-0.11).
VIF of dStromal_mid given ePV_2 + Arterial + uSMCs = 1.19 (R2 0.16): linear collinearity is low.
R2 pairwise conditional (p of row | Average + column):
- KOLLER_EUR: dStromal_mid 0.012 -> 0.078 | ePV_2; ePV_2 0.014 -> 0.093 | dStromal_mid (mutual
  attenuation; neither dominates); Arterial 0.007-0.017 unaffected by any other type; uSMCs 0.021
  -> 0.076 | dStromal_mid.
- GCST90205183: no type < 0.05 in any model (dStromal_mid 0.068-0.118).
- HYPERTENSION: Arterial, mPV, ePV_1b p < 0.001 robust to all conditioning; dStromal_mid 0.2-0.8.
- T2D: ePV_1b 0.002, mPV 0.007, Arterial 0.013, dStromal_mid 0.022 (-> 0.058 | ePV_2): a pattern
  similar to endometriosis.
- LEIOMYOMA: uSMCs p < 0.001 robust to every conditioning.
R3 delete-one-block jackknife (50 blocks): jackknife SEs close to MAGMA SEs (KOLLER dStromal_mid
0.46 vs 0.44; G1 lymphocyte 0.021 vs 0.022). Within trait BETA_STD(dStromal_mid) - BETA_STD(ePV_2):
KOLLER 0.0002 (z 0.03), GCST z 1.13, T2D z 0.66, HYPERTENSION z -2.21 (prefers ePV_2).
Cross-trait (block-paired): dStromal_mid KOLLER - T2D -0.04 (z -0.07), KOLLER - HTN z 0.92;
Arterial GCST - HTN z -2.86 (p 0.004), KOLLER - HTN z -2.00; G1 lymphocyte GCST - T2D z -1.91
(p 0.056); G1 epithelial no differences.
R4 simulations (null = KOLLER gene z circularly shifted; 30 replicates per scenario; observed
standardized effect 0.0196 per SD):
- null: dStromal_mid p < 0.05 in 3%, Holm-3 3%, BH-32 3% (calibrated).
- truth dStromal_mid 1x: detected p < 0.05 50%, Holm-3 27%, ranked first 13% (median rank 5);
  2x: Holm-3 97%, first 70%.
- truth ePV_2 1x: dStromal_mid Holm-3 7%, ePV_2 first 20%; 2x: dStromal_mid p < 0.05 53%, Holm-3
  37% (misassignment), never ranked first; ePV_2 first 70%.
- truth Arterial: dStromal_mid Holm-3 0-3% (separable). truth uSMCs 1x: dStromal_mid Holm-3 13%.
R5 true null (random patient splits within one context, KOLLER, 150 splits): P(p < 0.05) per
direction 0.19 / 0.11 (eutopic stromal, 9 patients), 0.25 / 0.21 (eutopic epithelial, 9), 0.12 /
0.13 (lesion stromal, 12); t SD 1.40, 1.81, 1.33 (expected 1).
B0: fallopian-tube atlases on CELLxGENE: Weigert et al. 2025 (131,218 cells, pre- and
post-menopausal, 1.6 GB); a post-menopausal tube + ovary atlas; Ulrich et al. 2022 (hydrosalpinx).
K0: no Koller supplementary file in references/.

Readings: (1) positive controls are resolved cleanly (leiomyoma -> uSMCs, hypertension -> arterial /
perivascular) while the endometriosis decidual signal cannot be separated from ePV_2 (z 0.03) and is
indistinguishable in size from T2D's; (2) at the observed effect size a true decidual source would
be detected (Holm-3) in only 27% of replicates and ranked first in 13% -> top-ranked cell types
from one GWAS are unreliable; a perivascular source at twice that size yields "significant"
decidual stroma in 37%; (3) gene-property tests of donor-estimated contrasts have 2-5x the nominal
false-positive rate even when both groups come from the same tissue (donor heterogeneity programs),
while specificity covariates are calibrated (R4 null 3%).

## 2026-11-07 · run p9_v028_scdrs_core_20261106T230906Z · Ochoa-faithful scDRS + donor-level tests

Cells (<= 400 per sample x label): 96,994 (peritoneal lesion 48,155; endometrioma 15,690;
peritoneum "no endometriosis detected" 10,465; eutopic 22,684). Runtime 975 min: later scoring runs
slowed ~60-fold (memory pressure across contexts on the laptop); results unaffected.
lesion_ochoa (74,310 cells; M2like = 4,305 of 8,611 myeloid cells):
- MC association, macrophage groups (myeloid, MPhi4, MPhi5, MACRO_ALL, M2like): GCST90205183
  p 0.073-0.134 (Holm 0.365); KOLLER_EUR 0.148-0.236 (Holm 0.739); HYPERTENSION 0.011-0.044
  (Holm 0.055-0.104); T2D 0.04-0.26; LEIOMYOMA 0.70-0.86.
- Donor-level z vs matched control sets, macrophage groups: GCST myeloid p 4e-07 (27/30 donors),
  M2like 1.5e-06 (23/26); KOLLER myeloid 3e-07 (29/30), M2like 5e-06 (23/26); HYPERTENSION myeloid
  3e-13 (30/30); T2D myeloid 1.4e-06 (27/30); LEIOMYOMA p 1.0 (2/30).
- Top labels: KOLLER eF2 p 0.001 (q 0.02), dS4 0.002 (q 0.02), eF3 0.002 (q 0.02); GCST eF3 0.004,
  eF2 0.006 (q 0.10); LEIOMYOMA eF2 0.001 (q 0.03); HYPERTENSION EC-capillary / EC-tip / pre-cDC2
  0.001 (q 0.01).
lesion_ours gives the same picture. Eutopic: KOLLER dS2 0.001, eF3 0.002, fib C7 0.002, eF2 0.003
(q 0.02); LEIOMYOMA eF3 0.003, fib C7 0.006; HYPERTENSION macrophages Holm 0.044.
Verdict by the fixed rule (v0.2.18): Ochoa claim NOT reproduced (MC Holm >= 0.05 for both
endometriosis GWAS). Reading: a donor-consistent excess of disease-gene expression in lesion
macrophages exists for endometriosis but equally for hypertension and T2D (not leiomyoma), i.e. it
is a generic property of complex-trait gene sets, not endometriosis-specific; the strongest
endometriosis scDRS signal is in endometrial-type stromal states (eF2, eF3, dS2, dS4), shared with
leiomyoma.

## 2026-11-07 · run p9_v030_koller_repro_20261107T160521Z · Koller 2026 method-matched reproduction

Cells: HECA dataset "Mareckova" (E-MTAB-14039), non-hormone: 64,791 cells; controls 7 donors
(26,327 cells), endometriosis 7 donors (38,464 cells). Koller report 22 endometriosis and 13 control
"samples" -> sample/donor correspondence not yet verified (v031 audit). <= 400 cells per type.
Markers: 91 sets (ABS and UP; all >= 20 genes); ABS median 1,829 genes. LD scores 20 min.
S-LDSC, family ABS (46 annotations incl. all-genes control), one-sided coefficient p, BH in family:
- KOLLER_EUR top: Venous/ENDO p 1.5e-04 (q 0.005), Venous/CTL 2.4e-04 (q 0.005), Arterial/CTL
  6.7e-04 (q 0.010), ePV_1a/CTL 1.4e-03 (q 0.016), eStromal/CTL 3.8e-03 (q 0.034),
  Glandular_secretory/ENDO 4.6e-03 (q 0.035).
  Focus: dStromal_early/ENDO p 0.017 (q 0.06), dStromal_mid/CTL 0.020 (q 0.06), dStromal_early/CTL
  0.040, SOX9_functionalis_II/CTL 0.042 (q 0.08), SOX9_functionalis_II/ENDO 0.072, dStromal_mid/ENDO
  0.31. (SOX9_functionalis_I and dStromal_late: no set in this cohort.)
- GCST90205183: Venous/ENDO 2.3e-03, Venous/CTL 2.5e-03 (q 0.055); dStromal_mid/CTL 0.009 (q 0.08).
- HYPERTENSION: Arterial/ENDO 0.0085 (q 0.38), nothing q < 0.05; T2D: nothing p < 0.09;
  LEIOMYOMA: nothing q < 0.05 (uSMCs not in top 6).
MAGMA gene-set (competitive, conditional on all genes):
- KOLLER_EUR ABS top Venous/CTL 3e-05, Arterial/CTL 9e-05, Venous/ENDO 1.3e-04; UP top
  dStromal_mid/ENDO 1.4e-04 (q 0.01), dStromal_early/ENDO 2.6e-04 (q 0.01).
- HYPERTENSION ABS Arterial/ENDO 6.7e-07, dStromal_early/CTL 2.8e-05; UP Arterial/ENDO 3e-08,
  dStromal_mid/ENDO 0.002. LEIOMYOMA UP uSMCs/CTL 1.5e-06 (positive control works in MAGMA).
Verdict by the fixed rule (v0.2.20): Koller nominal dStromal + SOX9 functionalis result REPRODUCED
(nominal); does NOT survive BH; "endometriosis-specific" YES by the rule (reproduced annotations
p >= 0.05 for HTN and T2D).
Reading: (1) the nominal claim reproduces, but under the identical pipeline the annotations that
survive correction are endothelial (venous, arterial), perivascular (ePV_1a), non-decidualized
stroma and secretory glandular, not decidual stroma; (2) the top cell type changes with the
selection rule: |log2FC| sets put endothelium first, up-only sets (MAGMA) put decidual stroma first
-> absolute-fold-change sets probably carry genes depleted in the small endothelial clusters, i.e.
genes of the majority stromal/epithelial cells (to be tested, v031); (3) the specificity reading is
weak: the S-LDSC positive controls failed (HTN -> arterial q 0.38; leiomyoma -> uSMCs absent) while
MAGMA detects them, so FinnGen comparators have little S-LDSC power (1000G EUR LD vs Finnish LD).
In MAGMA the decidual-stroma sets are also associated with hypertension.
B1 (v029) completed in the same session; its log is still to be read.

## 2026-11-07 · run p9_v029_b1_ft_20261107T160225Z · B1 fallopian tube (exploratory) — stopped at the compartment step

Weigert 2025 CELLxGENE atlas downloaded (1,646 MB; 131,218 cells, 17 donors, 10x 3' v3 + Drop-seq,
all "normal"); pre-menopausal by donor_menopausal_status: 85,107 cells; cell-type column
author_cell_type (70 labels); eligible pre-menopausal types 32 of 41; 13,747 genes shared with HECA.
Descriptive test (i), FT types | Average_FT (MAGMA gene-property, direction greater), top of list:
- KOLLER_EUR: st3 (stroma) p 0.010, st1 0.045, sm (smooth muscle) 0.053.
- GCST90205183: sm 0.034; nothing else < 0.1.
- ENDO_OVARY: st3 2.5e-04, sm 0.025, st2 0.026, ma 0.034, en2 0.040. ENDO_PERITONEUM: none < 0.05.
- HYPERTENSION: le (lymphatic endothelium) 7.7e-04, en1 (endothelium) 1.3e-03, p/v1 (perivascular)
  1.5e-03 -> vascular positive control works in this atlas.
- T2D: st3 5.6e-03, sm 0.025. LEIOMYOMA: tnk_pre5 0.007, sm 0.010.
No secretory (se_*) or ciliated (ce_*) type is among the top 5 for any endometriosis GWAS.
Compartment tests (ii)/(iii) NOT run: the pre-registered label patterns matched none of the author
labels (se_*, ce_*), the test table was empty and the script stopped (KeyError); the completion
marker was printed anyway (fixed in v032).
Reading (descriptive only): in the tube, as in the uterus, the endometriosis signal sits in stromal
programs; the same stromal type is also first for T2D.

## 2026-11-07 · run p9_v032_b1_ft_20261107T200330Z · B1 fallopian tube, re-run with amended labels (exploratory)

Label check (pre-menopausal cells): FT_SEC = se_pre2/3/4/5/7/8/10 (9,818 cells; CELLxGENE "peg cell"
83%, "fallopian tube secretory epithelial cell" 17%); FT_CIL = ce_pre2/3/5/7 (7,121 cells; 100%
ciliated epithelial cell); FT_STR (post hoc) = st1/st2/st3 (33,815 cells; 100% stromal cell).
HECA_EPI = 11 non-ciliated epithelial types, HECA_CIL = Ciliated, preCiliated, HECA_STR = 6 stromal
types. 13,747 shared genes. MAGMA gene-property | Average_FT + Average_HECA (one-sided p):
- KOLLER_EUR: FT_SEC 0.755 (| HECA_EPI 0.433); HECA_EPI 0.941; FT_CIL 0.784 (| HECA_CIL 0.25);
  HECA_CIL 0.937. Jackknife differences: SEC-EPI +0.0086 (se 0.0069, p 0.107); CIL-CIL +0.0066
  (se 0.0049, p 0.089).
- GCST90205183, ENDO_OVARY, ENDO_PERITONEUM: no epithelial compartment p < 0.05 marginally.
- Controls: no epithelial compartment p < 0.05 for HYPERTENSION, T2D or LEIOMYOMA.
Verdict (v0.2.19 rule): FT_SEC preferential over HECA_EPI: NO; FT_CIL over HECA_CIL: NO.
-> No evidence that endometriosis risk genes are preferentially specific to tubal secretory or
   ciliated epithelium; epithelium of either organ carries no detectable signal.
POST HOC stromal pair (hypothesis-generating): KOLLER_EUR FT_STR p 0.037 (| HECA_STR 0.155),
HECA_STR p 0.046 (| FT_STR 0.206), jackknife difference +0.0008 (se 0.0074, p 0.46) -> tubal and
uterine stromal programs are statistically indistinguishable and attenuate each other.
ENDO_OVARY: FT_STR 0.0057 (| HECA_STR 0.017), HECA_STR 0.087 (| FT_STR 0.53), difference p 0.094.
T2D: FT_STR 0.084, HECA_STR 0.089 (same direction, weaker). Verdict (same rule): NO.
Reading: the stromal association does not resolve to an organ; consistent with the HECA results
(decidual stroma not separable from perivascular, shared with T2D).

## 2026-11-07 · v031 interrupted (user, Ctrl-C) during P_UP LD scores; partial output

P markers re-derived: ABS 45, UP 45, DOWN 45 sets (median 1,829 / 644 / 981 genes); ABS sets
identical to v030 -> v030 LD scores reused. P_ABS adds ENDO_ALL (FinnGen R13 endometriosis; nested
in GCST90205183, probably also in KOLLER_EUR -> not independent): Venous/CTL p 6e-04 (q 0.03),
Venous/ENDO 2e-03 (q 0.04), Arterial/CTL 3e-03 (q 0.04); dStromal_early/ENDO 0.013,
dStromal_mid/CTL 0.019 -> a FinnGen-only GWAS with the same LD mismatch as the FinnGen control
traits does detect these annotations in S-LDSC. Full v031 readings pending (resumed by v033).

## 2026-11-07 · run p9_v031_koller_resolution_20261107T201721Z · Koller pipeline: selection rule, joint resolution, sensitivity (complete; resumed via v033)

K1 cohort (sample column "library"): Mareckova non-hormone ENDO 7 donors / 12 libraries, CTL 7 / 13;
with hormone-treated samples ENDO 9 / 14, CTL 7 / 13; all HECA non-hormone ENDO 19 / 24, CTL 30 / 41.
Koller: 22 endometriosis / 13 control samples -> the control count matches exactly, the
endometriosis count does not (12-14 vs 22; source of the difference unknown, possibly nuclei).
K3 S-LDSC families (P cells; ABS 45, UP 45, DOWN 45 sets; median 1,829 / 644 / 981 genes):
- KOLLER_EUR ABS: Venous/ENDO 2e-04 (q 0.01), Venous/CTL 2e-04 (q 0.01), Arterial/CTL 7e-04 (q 0.01).
  UP: dStromal_mid/ENDO 0.003 (q 0.12), dStromal_mid/CTL 0.013, dStromal_early/CTL 0.017; no
  endothelial set in top 6. DOWN: Ciliated/ENDO 6e-04 (q 0.03), Venous/ENDO 0.002 (q 0.04),
  Venous/CTL 0.002 (q 0.04).
- ENDO_ALL (FinnGen; nested): ABS Venous/CTL q 0.03; UP dStromal_mid/ENDO and /CTL 0.003 (q 0.05);
  DOWN Venous q 0.05.
- HYPERTENSION, T2D, LEIOMYOMA: nothing q < 0.19 in any family.
K2 composition: down share of ABS sets: endothelial Arterial/CTL 61%, Venous/CTL 57%, Venous/ENDO
61%, Arterial/ENDO 27%; all other sets median 60% (n = 41); decidual-stroma focus sets 82-88%,
SOX9_functionalis_II 53-74%. Jaccard focus x top sets median 0.25 (max 0.35).
H-ABS verdict: NOT SUPPORTED (c1 NO: Arterial/ENDO 27%; c2 YES; c3 YES). Down-dominance is not
specific to endothelium: most |log2FC| sets consist mainly of genes depleted in the labelled type;
the decidual-stroma sets that reproduce Koller's nominal result are 82-88% depleted genes.
K4 joint S-LDSC (focus | other / other | focus):
- ABS, KOLLER_EUR: 24 eligible pairs: "other only" 22, neither 2. No focus set survives
  conditioning on any top set (p 0.07-0.35); Venous/ENDO, Venous/CTL, Arterial/CTL, ePV_1a/CTL and
  Gland_sec/ENDO survive conditioning on every focus set. GCST90205183: other only 9, both 3,
  neither 2, focus only 1; no focus set survives all.
- UP, KOLLER_EUR (columns eStromal/CTL, sHormones/CTL, eHormones/CTL, eStromal_MMPs/CTL and /ENDO,
  Arterial/ENDO): 16 eligible: both 6, focus only 6, neither 4; dStromal_mid/CTL and /ENDO survive
  every column (p <= 0.04); only eHormones/CTL survives every focus set. GCST90205183: 4 eligible,
  all both; all four dStromal sets survive.
K6: ENDO_ALL reproduces the Koller sets (dStromal_early/ENDO 0.013, dStromal_mid/CTL 0.019,
dStromal_early/CTL 0.035); no S-LDSC positive control detected in any family (HTN best ePV_1b/ENDO
q 0.19; leiomyoma best uSMCs/CTL q 0.40-0.49) -> SPECIFICITY_VERDICT: comparators not informative.
K5 sensitivity (KOLLER_EUR, ABS): P nominal YES / BH NO; H (hormone-treated included; ENDO 9 donors)
nominal YES / BH YES (dStromal_late/ENDO 0.0067 q 0.03, dStromal_early/ENDO 0.016 q 0.04); N (all
cells) nominal YES / BH NO (SOX9f_II/CTL 0.0085, dStromal_mid/CTL 0.025). Top annotation in all
three variants: venous endothelium (q 0.01-0.03).
Reading: (1) the cell type "prioritized" under Koller's own pipeline is venous endothelium in every
variant; (2) with |log2FC| sets the decidual-stroma annotations are explained by the endothelial /
perivascular annotations and not vice versa; with up-only sets the direction reverses and
decidual stroma (mid) survives every conditioning -> the answer depends on the gene-selection rule;
(3) the |log2FC| annotations are mostly genes depleted in the named cell type; (4) specificity
cannot be judged with S-LDSC and the available comparators.

## 2026-11-07 · run p9_v034_koller_sign_20261107T224006Z · sign decomposition and decidual stroma vs vascular (v031 LD scores)

K7 (joint S-LDSC of each type's UP and DOWN sets; KOLLER_EUR), selected rows (ABS p | UP p, UP|DOWN | DOWN p, DOWN|UP):
- Venous/ENDO 2e-04 | 0.107, 0.050 | 0.002, 9e-04 -> DOWN-carried; Venous/CTL 2e-04 | 0.136, 0.069 |
  0.002, 0.001 -> DOWN-carried; Arterial/CTL 7e-04 | 0.153, 0.076 | 0.005, 0.003 -> DOWN-carried;
  ePV_1a/CTL 0.001 | 0.574 | 0.004, 0.003 -> DOWN-carried; mPV/CTL 0.049 -> DOWN-carried;
  Arterial/ENDO (27% down genes) 0.026 | 0.057, 0.046 | 0.132 -> UP-carried; ePV_2/CTL, uSMCs/CTL,
  SOX9f_II/CTL -> neither.
- dStromal_early/ENDO 0.017 | 0.032, 0.030 | 0.582 -> UP-carried; dStromal_mid/CTL 0.020 | 0.013,
  0.013 | 0.593 -> UP-carried; dStromal_early/CTL 0.040 | 0.017, 0.016 -> UP-carried;
  dStromal_mid/ENDO ABS 0.307 but UP 0.003 (UP|DOWN 0.003): signal masked in the |log2FC| set.
- ABS-nominal types: KOLLER_EUR 29 (DOWN-carried 15, UP-carried 7, neither 7); GCST90205183 19
  (DOWN 11, UP 3, both 1, neither 4); ENDO_ALL 23 (DOWN 15, UP 6, neither 2).
H-SIGN verdict (rule v0.2.23): NOT SUPPORTED (endothelial DOWN-carried 3/4 -- Arterial/ENDO is
UP-carried; decidual stroma UP-carried 3/3).
K8 (UP sets; dStromal x 11 vascular / perivascular / smooth-muscle sets, both directions): in all
three GWAS every dStromal set (early/mid x CTL/ENDO) survives conditioning on all 11 sets
(KOLLER p 0.003-0.032), and no vascular set survives conditioning on any dStromal set.
dStromal_mid | ePV_2: KOLLER 0.003-0.014, GCST 0.008-0.022, ENDO_ALL 0.004-0.027 -> G4 prediction
(not separable) NOT CONFIRMED in S-LDSC with up-only sets (G4 used MAGMA continuous specificity).
Reading: in Koller's |log2FC| layout more than half of the nominal annotations (15/29) get their
signal from genes depleted in the named cell type; the endothelial "top hits" are carried by
depleted genes; decidual-stroma signal is carried by up-regulated genes and, with up-only sets,
is separable from every vascular program (nominal level; marginal UP q 0.12 KOLLER, 0.05 ENDO_ALL).
Separability of decidual stroma from perivascular ePV_2 is framework-dependent (MAGMA specificity:
not separable; S-LDSC up-only markers: separable).

## 2026-11-07 · run p9_v035_abs_sim_20261107T231301Z · known-source simulation, design 1 (z circularly shifted)

19,261 genes; HECA profiles for 14,337; effect x1 = 0.0196 z per SD of profile. Observed KOLLER_EUR
(MAGMA gene-set): ABS top Venous/CTL, UP top dStromal_mid/ENDO, DOWN top Luminal/ENDO. 60 replicates
per scenario; shares of replicates with the top set in each group (dS / endothelial / pv+SM / other):
- null: ABS 0.08 / 0.07 / 0.17 / 0.68, any q<0.05 0.17; UP 0.15 / 0.10 / 0.12 / 0.63, any q 0.27;
  DOWN any q 0.17 -> NOT calibrated by the rule (<= 0.10).
- dStromal_mid x2: ABS 0.05 / 0.02 / 0.03 / 0.90 (median rank of best dS set 18); UP 0.27 / 0.00 /
  0.12 / 0.62 (median rank 3; dS q<0.05 0.28); x1: UP dS 0.17.
- Venous x2: ABS endothelial top 0.23 (median rank 4); UP 0.62 (rank 1; endothelial q<0.05 0.40).
- ePV_2 x2: ABS pv/SM top 0.13 (null 0.17); UP 0.63.
Verdict (v0.2.24): not calibrated; ABS relabels a decidual source: NO (endothelial/perivascular top
0.05 vs null 0.23); UP recovers the decidual source: NO (0.27 < 0.5).
Reading: under a planted source the |log2FC| rule rarely recovers it (5% / 23% / 13% for decidual /
venous / perivascular at 2x) while up-only sets recover vascular sources well (62% / 63%) and the
decidual source partly (27%); with a decidual source the ABS top set is "other" in 90%. The null is
anti-conservative in this design, most likely because shifting z against fixed gene-gene
correlations (stored in genes.raw, used by MAGMA) breaks the correlation structure; the comparison
between rules shares the same background. v036 repeats the simulation with the annotation layer
rotated instead (real z and correlations kept).

## 2026-11-08 · run p9_v036_abs_rot_20261108T000041Z · known-source simulation, design 2 (annotation layer rotated)

Real KOLLER_EUR z and gene-gene correlations; marker sets and planted profile rotated together;
60 replicates per scenario. Shares of replicates by top-set group (dS / endothelial / pv+SM / other):
- null: any q<0.05 ABS 0.02, UP 0.03, DOWN 0.02 -> CALIBRATED (v035 failure attributed to the
  shifted-z design, as pre-registered).
- dStromal_mid x2: ABS 0.02 / 0.02 / 0.12 / 0.85 (most often eHormones/ENDO, eStromal_MMPs/CTL,
  eStromal/CTL; best dS rank 19); UP 0.37 / 0 / 0.08 / 0.55 (eStromal/ENDO 0.30, dS_early/ENDO
  0.15; best dS rank 2; any q 0.28). x1: UP dS 0.23. x4: ABS dS 0.00, top preGlandular/ENDO 0.28,
  eStromal_MMPs/CTL 0.23 (dS q<0.05 0.00); UP any q 0.97, dS q<0.05 0.90, but top eStromal/ENDO
  0.40, dS_mid/ENDO 0.20, eStromal/CTL 0.17 (dS group 0.32); DOWN any q 0.75 (top Cycling,
  eHormones, Ciliated; endothelial q<0.05 0.18).
- Venous x2: ABS endothelial top 0.72 (Arterial/ENDO 0.37, Venous/CTL 0.17, Venous/ENDO 0.12);
  UP 0.95 (Venous/ENDO 0.52, Venous/CTL 0.30).
- ePV_2 x2: ABS pv/SM top 0.15 (null 0.17); UP 0.67 (ePV_2/CTL 0.30, ePV_2/ENDO 0.18).
Verdict (v0.2.25): calibrated YES; ABS relabels a decidual source: NO (0.13 vs null 0.27); UP
recovers the decidual source: NO (0.37 < 0.5).
Reading: (1) with correct correlation structure the competitive tests are calibrated (2-3%);
(2) |log2FC| sets recover a venous source (72%, mostly via the up-dominated Arterial/ENDO set)
but not a perivascular (15% vs null 17%) or decidual source (0-2%), which they attribute to
other stromal / glandular-progenitor states; (3) up-only sets recover vascular sources (95%, 67%);
a strong decidual source is DETECTED in 90-97% of replicates but ranked first only in 32-37%,
mostly losing to non-decidualized stroma (eStromal) -> detection without resolution among
stromal states; (4) under a decidual source, DOWN sets of other types (cycling, ciliated) become
significant (75% at x4) -> depleted-gene annotations can carry a stromal source.
Neither design reproduces the observed |log2FC| top hit (venous endothelium, DOWN-carried) from a
decidual source; v037 tests a lineage-agnostic parenchymal source.

## 2026-11-08 · run p9_v037_abs_par_20261108T005008Z · known-source simulation, parenchymal source (rotation design)

PAR profile = summed HECA specificity over 20 stromal + epithelial types (Ciliated, Glandular*,
HOXA13, KRT5, Luminal, MUC5B, SOX9*, dStromal*, eStromal*, preCiliated, preGlandular, preLuminal).
Null identical to v036 (calibrated: any q<0.05 ABS 0.02, UP 0.03, DOWN 0.02).
- PAR x2: ABS top dS 0.15 / endothelial 0.03 / pv+SM 0.03 / other 0.78 (preCiliated/ENDO,
  Ciliated/CTL); UP other 0.93 (Ciliated/CTL 0.22); DOWN endothelial 0.23, pv+SM 0.28
  (Venous/ENDO 0.10, ePV_2/ENDO 0.08).
- PAR x4: ABS endothelial 0.00 (top preCiliated/ENDO 0.22, Ciliated/CTL 0.20); UP Ciliated/CTL
  0.50; DOWN endothelial 0.32, pv+SM 0.45 (endothelial q<0.05 0.40).
- PAR x2 + dStromal_mid x2: ABS other 0.92; UP dS 0.38; DOWN endothelial 0.25.
H-PAR verdict (v0.2.26): NOT SUPPORTED (ABS top endothelial: null 0.10, PAR x2 0.03, x4 0.00).
Observed KOLLER signature (ABS top endothelial AND UP top stromal): 0.00 in every planted scenario,
0.03 under the null.
Reading: a parenchymal source lights up the DOWN sets of vascular types (as expected) but not the
|log2FC| vascular sets, which are outranked by ciliated sets; the summed-specificity profile is
dominated by sharply specific small types (ciliated), so it is not the lineage-agnostic source
intended -- a limitation of the H-PAR test. No single planted source tested (decidual,
venous, perivascular, parenchymal, mixture) reproduces the observed |log2FC| venous top hit; the
real-data mechanism stands as measured (K7: the venous |log2FC| signal is carried by genes
depleted in venous endothelium), its source configuration unresolved.

## 2026-11-08 · v038 draft figures (no analysis)

F1 rule flip, F2 sign decomposition, F3 joint resolution, F4 simulations, F5 |log2FC| composition
written to results/paper_figures/p9_v038_*/ (PNG, PDF, source-data CSV).

## 2026-11-08 · run p9_v039_locus_atlas_20261108T014217Z · locus atlas (partly invalid: gene matching bug)

Loci (distance clumping, p < 5e-8): KOLLER_EUR 57, GCST90205183 24, leiomyoma 138, hypertension 408,
T2D 375. BUG: for KOLLER_EUR and the FinnGen traits no candidate gene was matched ("no expression"
for every locus) because the genes.out chromosome column was read as text (chrX rows); only
GCST90205183 (no chrX genes) was processed: 24 loci, 13 assigned -- endothelium 4, immune 4, stroma
2, perivascular/SM 2, epithelium 1 (9 broad, 2 without expressed genes). L1 and L3 compartment
results for the other traits are void; v040 re-runs them.
Valid (independent of the bug): L2 leave-one-locus-out on KOLLER_EUR (57 loci; gene removal by
gene.loc positions): baseline ABS top Venous/CTL 3.2e-05, Arterial/CTL 8.7e-05, Venous/ENDO 1.3e-04;
UP top dStromal_mid/ENDO 1.4e-04, dStromal_early/ENDO 2.6e-04, sHormones/CTL 9.0e-04; DOWN top
Luminal/ENDO 5.3e-05, Ciliated/CTL 4.8e-04, Ciliated/ENDO 5.1e-04. No single locus raises the p of
any tracked set by >= 10x -> the aggregate rankings are not driven by individual loci.
L3 counts (45/57 loci "shared" with at least one of leiomyoma / T2D / hypertension at p < 1e-5)
lack a baseline (hypertension and T2D have tens of thousands of significant SNPs); v040 adds
random windows.

## 2026-11-08 · run p9_v040_locus_atlas_20261108T015742Z · locus atlas (corrected gene matching)

K1b: E-MTAB-14039 sample table: 103 rows, 63 libraries, 74 individual ids; disease endometriosis 71
rows, normal 32; stage Hormones 18 rows; all rows "10x 3' v1, single cell". All 22 HECA Mareckova
(cells file) libraries are in the table; they cover C 10 libraries (6 donors), E 12 libraries
(hormones 2, E + F 2). Libraries pool several donors (library ids repeat across individuals), so
"samples" is ambiguous; K1c (v042) tabulates individuals / pairs / libraries.
L1 (57 KOLLER_EUR loci; 52 with an expressed candidate gene): assigned 32 -- stroma 9 (28%; WNT4,
GREB1, PARPBP/IGF1, TBX3, WT1, TSHZ3, HOXA9-11, VPS13B/OSR2, SFRP1), endothelium 8 (25%; KDR, FGD6/VEZT,
CD109, RGS17, SRP14, DNM3, CALCRL, SKIDA1), immune 6 (19%; PTPRO, SKAP1, RIN3, OPRM1, CKAP2L/IL1A/IL1B,
HLA-DRA region), epithelium 5 (16%; ABO, HMGA1, WNT7A, BMPR1B, RUNX1), perivascular/SM 4 (SYNE1 x2,
RSPO3, NGF); broad 20 (incl. FSHB, ESR1, MTAP/CDKN2A-B, ID4, GNRH1); no expressed gene 5.
Enrichment vs all expressed genes: stroma 9/52 (genome 8%) p 0.030; endothelium 8/52 p 0.069; others
n.s. H-MULTI: SUPPORTED (5 compartments >= 10%, top share 0.28) but not distinctive: leiomyoma (0.25),
hypertension (0.24) and T2D (0.23) are equally spread.
L2: no single locus raises any tracked set's p by >= 10x. L2b (all 57 loci removed, 468 genes):
ABS Venous/CTL 3.2e-05 -> 0.14, Arterial/CTL 8.7e-05 -> 0.18, Venous/ENDO 1.3e-04 -> 0.17; UP
dStromal_mid/ENDO 1.4e-04 -> 0.021, dStromal_early/ENDO 2.6e-04 -> 0.086; DOWN Ciliated/CTL 4.8e-04 ->
0.017. Verdict: 3/9 tracked sets keep p < 0.05 -> "signal mostly in genome-wide loci" (jointly).
New tops without the loci: ABS eStromal_MMPs 0.057; UP dStromal_mid/ENDO 0.021.
L3 (unmatched baseline, 382 random windows): shared with leiomyoma 26/57 (46%) vs 8% (OR 9.8, p 1e-11);
T2D 24/57 (42%) vs 27% (OR 2.0, p 0.016); hypertension 28/57 (49%) vs 29% (OR 2.4, p 0.002).
Matched baseline in v042.
Reading: the |log2FC| endothelial hit depends jointly on the genome-wide loci (no single one), while
the up-only decidual-stroma signal partly persists without them; several of the strongest loci act
on the endocrine axis outside the endometrium (FSHB, GNRH1; ESR1 broad), which no endometrial atlas
can localize.

## 2026-11-08 · run p9_v041_benchmark_20261108T020555Z · matched directional benchmark (conditional; rotation design)

Thresholded sets identical to v031 (136/136). 270 sets (6 families x 45) + 45 continuous covariates.
Calibration (300 null rotations; FWER = any BH q < 0.05 within family; 95% Clopper-Pearson):
ABS 0.010 (0.002-0.029), UP 0.003, DOWN 0.007, ABSN 0.003, UPN 0.007, DOWNN 0.003, SPEC 0.017
(0.005-0.038) -> all consistent with <= 5% (conservative).
Top-1 exact type recovery (40 replicates per truth; strength 2x):
- POS (specificity): Venous SPEC 0.93, UP 0.65, UPN 0.60, ABSN 0.55, ABS 0.12, DOWN 0.00;
  dStromal_mid SPEC 0.47 (top-3 0.80, group 0.90), UP 0.20, UPN 0.15, ABSN 0.03, ABS 0.03;
  ePV_2 SPEC 0.78, UP 0.50, UPN 0.47, ABSN 0.15, ABS 0.12.
- NEG (depletion): only down-sets recover: Venous DOWN 0.28 / DOWNN 0.15; dStromal_mid DOWN 0.20 /
  ABS 0.07; ePV_2 DOWN 0.15 / DOWNN 0.07; UP, UPN, SPEC 0.00 everywhere.
- EXPR (expression level): best SPEC (Venous 0.42, ePV_2 0.23, dStromal_mid 0.15); ABS 0.00.
- Related-state misses for the decidual truth: UP / UPN correct group 0.82 but exact 0.15-0.20
  (miss 0.62-0.68); SPEC group 0.90, exact 0.47 (miss 0.42).
Verdicts (v0.2.30): P1 YES (UPN >= ABSN under POS for all targets); P2 YES (down / |log2FC| beats
up under NEG); P3 NO (UPN - ABSN >= 0.2 only for ePV_2: 0.32; Venous 0.05, dStromal_mid 0.12);
P4 (exploratory) YES (no family > 0.5 under EXPR).
Reading: no rule is best everywhere; the thresholded |log2FC| rule (Koller) is never the best and
is the worst or near-worst in every scenario; much of its deficit is set size / composition
(ABS 0.12 -> ABSN 0.55 for Venous); the direction of the true association decides which
directional rule works (up-type sources: SPEC / UP; depletion-type: DOWN only); decidual stroma is
located to the stromal group but not to the exact state even by the best method. Conditional
benchmark: truths built from expression features (SPEC advantage under POS partly by design).

## 2026-11-08 · run p9_v042_locus_atlas_20261108T042051Z · locus atlas, final form

K1c (Koller cohort): E-MTAB-14039 has 103 (library, individual) rows, 74 individuals, 63 libraries.
"Libraries in the HECA cells file, hormone-treated included" gives (library, individual) pairs
E 22 / C 13 = Koller's 22 endometriosis / 13 control "samples" EXACTLY (individuals 9 / 7; libraries
14 / 13). Hormone-excluded (our v030 primary cells): pairs 18 / 13. -> Koller's cohort = our v031
variant H (hormone-treated included): there the nominal decidual result also survives BH
(dStromal_late/ENDO q 0.03, dStromal_early/ENDO q 0.04) and the top annotation is still venous
endothelium (q 0.01).
L1 checks: KOLLER_EUR genes.out 19,261 genes (787 chrX excluded); 57 loci, 53 with candidate genes
(258 genes), 52 with an expressed gene; leiomyoma 129 loci (1,231 chrX SNPs excluded), hypertension
393, T2D 349, GCST 24. Compartments as in v040.
L2c (loci removed by compartment; KOLLER_EUR; MAGMA, v031 P sets):
- stromal loci (9: WNT4, GREB1, PARPBP/IGF1, TBX3, WT1, TSHZ3, HOXA9-11, VPS13B/OSR2, SFRP1): UP
  dStromal_mid/ENDO 1.4e-04 -> 5.2e-03 (+1.57 log10) vs random 9-locus subsets median +0.08, 95th
  +0.51 -> permutation p 0.020 -> SUPPORTED. The same removal also weakens ABS Venous/CTL
  (3.2e-05 -> 7.1e-04, +1.35), ABS Venous/ENDO (+0.79), DOWN Ciliated/ENDO (+1.37), DOWN
  Luminal/ENDO (+1.69).
- endothelial loci (8: KDR, FGD6/VEZT, CD109, RGS17, SRP14, DNM3, CALCRL, SKIDA1): ABS Venous/CTL
  +0.89 vs random 8-locus subsets median +0.68, 95th +1.21 -> permutation p 0.275 -> NOT SUPPORTED.
- broad loci (20, incl. FSHB, ESR1): ABS Venous/CTL -> 1.4e-03 (+1.64); UP dStromal_mid -> 3.5e-04.
L2b wording: "signal depends on the genome-wide loci jointly" (3/9 tracked sets keep p < 0.05).
L3 vs matched windows (pool 800, 399 controls; balance: median 3 genes / ~4,080 SNPs in loci,
matched and pool): leiomyoma 26/57 (46%) vs 39/399 (10%), OR 7.7, p 4e-10; T2D 42% vs 28%, OR 1.9,
p 0.022; hypertension 49% vs 30%, OR 2.3, p 0.0035.
Reading: the up-only decidual-stroma signal is carried specifically by a hormone-responsive
developmental stromal program (WNT4, GREB1, IGF1, WT1, HOXA10/11, TBX3, SFRP1); the |log2FC| venous
"top hit" is not carried by the endothelial-gene loci and is reduced at least as much by removing
the stromal loci (formal test in v044) -- consistent with the depleted-gene mechanism (K7).
Endometriosis loci are strongly shared with leiomyoma.

## 2026-11-08 · run p9_v043_real_families_20261108T043912Z · seven families on real GWAS, positive controls

E1 top per family (p / BH q within family):
- KOLLER_EUR: ABS Venous/CTL 3e-05; UP dStromal_mid/ENDO 1.4e-04 (q 0.01); DOWN Luminal/ENDO 5e-05;
  ABSN eStromal_MMPs/ENDO 7e-05; UPN dStromal_mid/ENDO 1.2e-06 (q < 0.01), dStromal_early/ENDO
  2.7e-06; DOWNN SOX9f_II/ENDO 3.6e-08; SPEC uSMCs/CTL 4.4e-03 (q 0.20), dStromal_mid/ENDO 0.020.
- GCST90205183: ABS Venous/CTL 1.2e-03 (q 0.03); UPN Arterial/CTL 8.9e-04 (q 0.04); SPEC
  dStromal_mid/ENDO 8.2e-03 (q 0.37).
- ENDO_ALL: ABSN ePV_2/ENDO 1.6e-06; UPN eStromal/CTL 6.5e-06; SPEC ePV_2/ENDO 0.010 (q 0.47).
- LEIOMYOMA: UP / UPN / ABSN / SPEC uSMCs/CTL first (UPN 8e-08, SPEC 5.7e-06); ABS dStromal_early.
- HYPERTENSION: ABS Arterial/ENDO 6.7e-07; UP Arterial/ENDO 3e-08; UPN Arterial/CTL 6e-06; SPEC
  mPV/CTL 8e-07; DOWN dStromal_early.
- T2D: SPEC dStromal_mid/ENDO 9.7e-05 (q < 0.01); UPN eStromal/CTL 1.5e-03 (q 0.07).
E2 verdicts (v0.2.33): leiomyoma -> smooth muscle first in UP, ABSN, UPN, SPEC (not ABS, DOWN, DOWNN);
hypertension -> vascular first in ABS, UP, ABSN, UPN, SPEC (not DOWN, DOWNN).
E3 (exploratory): no v041 scenario reproduces the full observed KOLLER pattern (0.00 jointly);
POS_dStromal_mid matches UP / UPN / DOWN / DOWNN in 0.82-0.90 of replicates but ABS only 0.03.
Reading: the thresholded |log2FC| rule fails the leiomyoma positive control that every up-type rule
passes; size-matched up-only gives the strongest endometriosis result (decidual stroma, 1e-06);
continuous specificity finds decidual stroma significant for T2D (q < 0.01) but not for KOLLER_EUR
(q 0.46), so "decidual stroma" by specificity is not endometriosis-specific.

## 2026-11-08 · run p9_v044_external_20261108T050151Z · stromal loci vs the venous hit; external atlas

A (H-STRVEN, KOLLER_EUR, v031 P sets; 9 stromal / 8 endothelial loci; 200 random 9-locus subsets):
- ABS Venous/CTL (primary): base 3.2e-05; without stromal loci 7.1e-04 (+1.35), without endothelial
  loci 2.4e-04 (+0.89); random median +0.69, 95th +1.19; perm p 0.020 -> SUPPORTED.
- ABS Arterial/CTL +1.17 vs +0.80 (perm p 0.070); ABS Venous/ENDO +0.81 vs +0.86 (perm p 0.318).
- DOWN Ciliated/ENDO +1.35 vs +0.15 (perm p 0.005); DOWN Luminal/ENDO +1.69 vs +0.28 (perm p 0.005);
  UP dStromal_mid/ENDO +1.57 vs +0.48 (perm p 0.005; 0.005 = minimum with 200 subsets).
Reading: the same nine stromal loci (WNT4, GREB1, IGF1 region, TBX3, WT1, TSHZ3, HOXA9-11, OSR2 region,
SFRP1) carry the up-only decidual-stroma signal, the epithelial down-set signals (luminal, ciliated)
and the primary |log2FC| venous hit: one stromal genetic signal appears under three cell labels
depending on the marker rule. Secondary vascular |log2FC| sets only partly follow (Arterial
directionally, Venous/ENDO not).
B (GTEx v9 snRNA-seq, Eraslan 2022; CELLxGENE dataset 4ed927e9-c099-49af-b8ce-a2652d069333; 209,126
nuclei; 8 tissues: heart LV, lung, oesophagus muscularis, prostate, gastrocnemius, oesophagus mucosa,
breast, skin): 45 of 53 cell types kept; 17,054 genes mapped. Wording fix: the atlas has no uterine
or female reproductive-tract tissue but includes breast and prostate, so "non-reproductive" is
inaccurate -> "an atlas without uterine or female reproductive-tract tissue".
- Positive controls (pre-registered): LEIOMYOMA -> smooth muscle first with UP, UPN, SPEC (ABS rank 7,
  ABSN 4, DOWN 16, DOWNN 8). HYPERTENSION -> SPEC only (cardiac endothelial first); in ABS, UP, ABSN,
  UPN the first type is "contractile cell" (Cell Ontology term; identity checked post hoc in v045) and
  cardiac endothelial is second. T2D (exploratory) -> subcutaneous adipocyte first with UP, UPN
  (ABS rank 10).
- Down-type families rank the same cell types first for every trait (squamous epithelial cell,
  cardiac fibroblast, slow/fast muscle; q < 0.01 for all six GWAS); ABS top hits coincide with DOWN
  top hits for KOLLER_EUR, LEIOMYOMA and T2D (cardiac fibroblast, fast muscle) -> formal tests v045.
- Endometriosis in this atlas: UP / UPN top = hormone-responsive epithelia (mammary basal-myoepithelial
  and luminal, prostate basal / luminal epithelium) and skin fibroblast (UPN prostate basal 8.4e-06);
  SPEC nothing at q < 0.10 (best q 0.14); T2D UPN also ranks mammary luminal second -> exploratory only.

## 2026-11-08 · run p9_v045_specificity_20261108T054316Z · label check, trait-specificity, gene-property conditioning

C1 (post hoc label check, GTEx v9): "contractile cell" (5,078 nuclei; heart 48%, oesophagus
muscularis 19%) = author Broad cell type Pericyte/SMC 98%, Granular Pericyte/SMC I/II 100%; marker
panels pericyte 938 / vascular SMC 1,001 / myoepithelial 2 mean CPM; top up genes PRKG1, DLC1, EBF1,
PDGFRB, CACNA1C. The fixed rule printed "not vascular mural" because it also required a mural
majority in the coarse hierarchy columns ("Cell types level 2" = Muscle, "level 3" = Stromal), which
never use mural vocabulary -> rule mis-specified; deviation declared in CHANGELOG v0.2.36.
"basal-myoepithelial cell of mammary gland" = 71% oesophageal-mucosa myofibroblasts + 29% breast basal
myoepithelium (same CL term) -> the v044 reading "hormone-responsive mammary epithelium" is withdrawn.
"smooth muscle cell" = oesophageal smooth muscle 95%; "enteric smooth muscle cell" = oesophageal
smooth muscle (TAGLN lo) 100%.
C2 H-NONSPEC (KOLLER_EUR, LEIOMYOMA, HYPERTENSION, T2D): mean cross-trait Spearman rho
- HECA: ABS 0.32, UP 0.43, DOWN 0.42, ABSN 0.32, UPN 0.47, DOWNN 0.34, SPEC 0.51; down-type minus
  up-type -0.09 (95% CI -0.28 to 0.11);
- GTEx: ABS 0.48, UP 0.38, DOWN 0.49, ABSN 0.47, UPN 0.36, DOWNN 0.58, SPEC 0.29; +0.20 (-0.00 to 0.39)
-> NOT SUPPORTED. Top-hit level (descriptive): HECA KOLLER ABS Venous/CTL ranks 10 / 4 / 26 in
leiomyoma / hypertension / T2D; UP dStromal_mid/ENDO 3 / 7 / 9; SPEC KOLLER top uSMCs/CTL ranks 1 for
leiomyoma and 2 for T2D, SPEC dStromal_mid/ENDO ranks 1 for T2D. GTEx KOLLER DOWN top (squamous
epithelial cell) ranks 1 for all three comparators; ABS top (cardiac fibroblast) 1 / 3 / 2.
C3 H-ABSDOWN: depleted-gene share of |log2FC| sets HECA median 0.60 (31/45 > 0.5), GTEx 0.33 (7/45);
rho(ABS, DOWN) < rho(ABS, UP) in 6/6 traits (HECA) and 5/6 (GTEx) -> NOT SUPPORTED. |log2FC| rankings
globally track up-set rankings more than down-set rankings; the depleted-gene mechanism (K7) applies
to specific hits, not to the whole ranking.
C4 (MAGMA gene-property analysis of binary set indicators; gnomAD v2.1.1 LOEUF 19,155 genes; HECA
14,205 genes, 270 indicators; GTEx 15,243 genes, 269 indicators): equivalence check with set analysis
Spearman 0.968 (HECA) / 0.984 (GTEx).
- BH-significant labels summed over 6 traits, base -> conditioned on Average + Breadth + LOEUF:
  HECA ABS 43->14, UP 16->7, DOWN 49->5, ABSN 34->8, UPN 32->10, DOWNN 69->7, SPEC 9->7 (retained
  down-type 0.10, up-type 0.42); GTEx ABS 143->20, UP 28->9, DOWN 127->49, ABSN 125->31, UPN 51->11,
  DOWNN 170->72, SPEC 4->1 (down-type 0.41, up-type 0.25) -> H-GENERIC NOT SUPPORTED (opposite in GTEx).
- Positive controls unchanged by conditioning: HECA leiomyoma first with UP, ABSN, UPN, SPEC (ABS
  2 -> 1); hypertension first with ABS, UP, ABSN, UPN, SPEC; GTEx leiomyoma first with UP, UPN, SPEC
  (ABS 7 -> 7).
- KOLLER_EUR key sets (HECA), base -> expression level only -> all covariates: ABS Venous/CTL
  1.2e-04 -> 1.2e-03 -> 7.6e-04 (rank 1 -> 1); ABS Arterial/CTL 2.3e-04 -> 1.6e-03; UP dStromal_mid/ENDO
  1.7e-04 -> 8.4e-04 -> 1.2e-03 (rank 1); UPN dStromal_mid/ENDO 1.5e-05 -> 5.3e-04 -> 8.9e-04 (rank 1);
  DOWN Luminal/ENDO 2.0e-05 -> 4.8e-04 (rank 1); SPEC dStromal_mid/ENDO 0.013 -> 0.030 (rank 2)
  -> H-DSTROMA SUPPORTED.
Reading: generic gene properties (expression level, breadth, LoF constraint) remove a large part of
the significance in both atlases but do not change which label each rule ranks first: the
rule-dependence of the top label (venous / decidual stroma / luminal) is not a gene-property artefact.
The "down-type families are trait-nonspecific" and "|log2FC| = down-set" generalisations are not
supported and are dropped.

## 2026-11-08 · run p9_v046_sign_20261108T151714Z · composition and sign decomposition of |log2FC| hits (MAGMA)

D1 Spearman(depleted-gene share, -log10 p of the |log2FC| set) across 45 sets: HECA +0.03 (KOLLER_EUR)
to +0.27 (hypertension); GTEx +0.29 to +0.51 (positive in all six traits).
D2 (K7 rule, joint MAGMA gene-property analysis of the UP and DOWN parts), top-5 |log2FC| hits with
p < 0.05 pooled over six traits: GTEx n 30 -> DOWN-carried 10, UP-carried 3, both 17 (share 0.33);
HECA n 27 -> DOWN-carried 10, UP-carried 9, both 7, neither 1 (0.37) -> H-K7REP NOT SUPPORTED.
Descriptive: the depleted part is independently significant (DOWN | UP < 0.05) in 27/30 (GTEx) and
17/27 (HECA) top hits; hits carried by the up part alone 3/30 and 9/27.
KOLLER_EUR key sets (HECA; ABS p | UP|DOWN | DOWN|UP): Venous/CTL 1.2e-04 | 0.012 | 9.6e-04 -> both
(S-LDSC K7: DOWN-carried); Arterial/CTL 2.3e-04 | 0.018 | 1.4e-03 -> both; Venous/ENDO 2.8e-04 |
5.5e-03 | 5.1e-03 -> both; Luminal/ENDO 2.5e-04 | 0.35 | 1.9e-05 -> DOWN-carried; dStromal_early/CTL
2.5e-03 | 0.10 | 5.2e-03 -> DOWN-carried (S-LDSC K7: UP-carried); Ciliated/ENDO -> DOWN-carried.
Masking (replicates the S-LDSC K7 observation): dStromal_mid/ENDO |log2FC| set p 0.35 (rank 41 of
45; 83% depleted genes) while its up part has p 1.7e-04 (| DOWN 2.2e-04); dStromal_early/ENDO 0.10 vs
up part 5.1e-04.
D3 (leiomyoma): HECA uSMCs/CTL |log2FC| p 0.027 (rank 2; 60% depleted) vs up part 4.7e-06; the one set
above it (dStromal_early/CTL, 88% depleted) DOWN-carried -> YES. GTEx "smooth muscle cell" |log2FC|
p 1.0e-03 (rank 7; 31% depleted; up part 0.021, down part 7.3e-03); six sets above: DOWN-carried 2,
both 4 -> NO (descriptive: in all six DOWN | UP < UP | DOWN, 3.9e-06 to 4.5e-03 vs 0.011 to 0.50)
-> H-MASK NOT SUPPORTED. Hypertension: GTEx Pericyte/SMC ("contractile cell") first under |log2FC|
(up part 6.9e-06); HECA Arterial/ENDO first (27% depleted; up part 4.6e-07).
Reading: masking of the true source by its own depleted genes is robust across frameworks (S-LDSC,
MAGMA) and atlases (uterine smooth muscle for leiomyoma: up part 4.7e-06 -> |log2FC| set 0.027;
decidual stroma for endometriosis: 1.7e-04 -> 0.35). Promotion of other cell types through depleted
genes is partial: the depleted part contributes in most top hits but usually together with the up
part; the pre-registered "mostly DOWN-carried" criterion is not met. Sign classes of individual sets
can differ between S-LDSC and MAGMA (Venous/CTL, dStromal_early/CTL).

## 2026-11-08 · run p9_v047_koller_exact_20261108T154620Z · exact Koller cohort (variant H) through S-LDSC

Sets: 49 |log2FC|, 49 up-only, 49 down-only (v031 variant H cells; hormone-treated samples included).
X1 / X2 first-ranked annotation (p / BH q within rule):
- KOLLER_EUR: |log2FC| Venous/ENDO 1e-04 / 0.01 (then Venous/CTL, Arterial/CTL); up-only
  dStromal_mid/ENDO 1.6e-03 / 0.08 (then sHormones/ENDO, dHormones/ENDO); down-only Venous/ENDO
  9e-04 / 0.04. Koller's decidual result under |log2FC|: dStromal_late/ENDO 0.0067 / 0.03.
- GCST90205183: |log2FC| Venous/ENDO 2e-03 / 0.05; up-only dStromal_mid/ENDO 6.5e-03 / 0.14.
- ENDO_ALL: |log2FC| Venous/CTL 6e-04 / 0.02; up-only sHormones/ENDO 8e-04 / 0.02, dStromal_mid/ENDO
  9e-04 / 0.02.
- LEIOMYOMA: uterine smooth muscle rank 8 with |log2FC|, 1 with up-only (uSMCs/CTL p 0.011, q 0.53),
  32 with down-only. HYPERTENSION: vascular rank 1 with |log2FC| (Arterial/ENDO p 0.021, q 0.43) and
  up-only (ePV_1b/ENDO 2.8e-03, q 0.14), 8 with down-only.
-> H-EXACT SUPPORTED (rank-based). The S-LDSC positive-control signals are nominal only (q >= 0.14).
K2: down share of |log2FC| sets: Venous/ENDO 62%, Venous/CTL 57%, Arterial/CTL 61%, Arterial/ENDO 36%;
other sets median 57%; decidual sets 59-88% -> H-ABS NOT SUPPORTED on variant H (Arterial/ENDO < 50%;
P: supported).
K4 |log2FC|: no decidual / SOX9-functionalis set survives any of the six top non-focus sets (KOLLER
33/36 pairs "other only"; GCST 17/20); every top non-focus set survives every focus set.
K4 up-only: dStromal_mid/ENDO and dStromal_mid/CTL survive all six top up-only sets (KOLLER);
dStromal_early/ENDO, dStromal_mid/CTL, dStromal_mid/ENDO (GCST).
K6: Koller-reproduced sets nominal in FinnGen ENDO_ALL (dStromal_early/ENDO p 0.004); no S-LDSC
positive control at q < 0.05.
K7 (KOLLER_EUR): 34 nominal |log2FC| annotations: DOWN-carried 19, UP-carried 8, both 1, neither 6
(variant P: 15/29 DOWN-carried). Venous/ENDO, Venous/CTL, Arterial/CTL, ePV_1a/CTL DOWN-carried;
Arterial/ENDO both; dStromal_early/ENDO, dStromal_mid/CTL, dStromal_early/CTL UP-carried;
dStromal_late/ENDO DOWN-carried. Masking: dStromal_mid/ENDO |log2FC| p 0.073 vs up part 0.002.
H-SIGN (all-of rule) NOT SUPPORTED (3/4 endothelial DOWN-carried; 3/4 decidual UP-carried).
GCST: 23 nominal: DOWN-carried 12, UP-carried 6, both 1, neither 4.
K8 (up-only): dStromal_mid/ENDO, dStromal_mid/CTL, dStromal_early/CTL, dStromal_early/ENDO survive
all 11 vascular / perivascular / smooth-muscle up-only sets (KOLLER); all five decidual sets survive
all 11 (GCST); dStromal_mid separable from ePV_2.
Reading: every variant-P result holds in Koller's exact cohort: under their rule the first-ranked
cell type is venous endothelium and the decidual-stroma signal is explained by endothelial /
perivascular annotations; under up-only selection decidual stroma ranks first and survives every
vascular program; the endothelial |log2FC| hits are carried by depleted genes; and the |log2FC| rule
ranks uterine smooth muscle 8th for leiomyoma, where up-only ranks it first (nominal in S-LDSC).
