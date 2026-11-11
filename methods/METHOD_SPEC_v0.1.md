# P9-DARC — Method Specification v0.1

**Working name:** Direction-aware and Resolution-calibrated Cell Prioritization (DARC).  
**Status:** Research design / NOT a validated method.  
**Date:** 2026-11-11  
**Relation to Paper 9:** A new *prospectively specified follow-up*, motivated by exploratory Paper 9 v001–v071. This is **not** retroactive preregistration.  
**Scope:** Inference on relative GWAS-to-cell-state associations, not proof of a causal cellular origin.

## 1. Scientific objective and target estimands

Existing methods such as CELLECT, scDRS and CATCH prioritize disease-relevant cell types. Our proposed **incremental claim to be tested** is that explicit marker *direction*, *donor uncertainty*, and *hierarchical resolution decisions* can reduce overconfident specific-cell assignments at usable coverage. Novelty is a hypothesis; do not claim superiority before external evaluation.

For each candidate cell state `c`, estimands are:

1. `beta(c, a, m)`: genome-wide enrichment of trait association in marker annotation `a` for cell `c`, using a designated genetics engine `m` (MAGMA gene-property or S-LDSC). Values from different engines are NOT interchangeable without standardized modeling.
2. `delta(c) = beta(c, UP, m) - beta(c, DOWN, m)`: a *within-engine, jointly estimated* directional contrast when both coefficients are identifiable and standardized. This is NOT a gene-expression effect of a risk allele.
3. `R_state(c)` and `R_group(parent(c))`: out-of-sample assignment-error risk among resolved state and group predictions, estimated in independently designed validation sets. They are NOT individual posterior probabilities of causality.
4. `coverage`: the fraction of tasks assigned a state or group rather than abstained.

Primary method success is **lower incorrect-state assignment risk at comparable coverage** than prespecified baselines, with error risk estimated under external, known-truth simulation regimes. A failed accuracy/coverage gate is informative and cannot be rescued by choosing a new threshold on the test set.

## 2. Input contract

### Required for the eventual complete method

- GWAS summary statistics (`chr`, `pos`, `effect allele`, `other allele`, effect/SE or z, sample size or effective N) with ancestry and reference build.
- An ancestry-appropriate LD reference, matching genome build and alleles.
- scRNA-seq donor-level data or a donor-by-cell-state pseudobulk matrix, plus donor IDs, tissue/context, experimental condition, batch, phase, and cell-state labels where available.
- A cell-state hierarchy (`state -> lineage/group`), prespecified from the atlas without examining the GWAS results.
- Gene identifiers mapped to one frozen genomic reference and one COMMON tested-gene universe for all annotation families.

### Optional

- Multiple independent reference atlases; independent ATAC/QTL data; scDRS/CELLECT/CATCH baseline results.
- If donor identities are unavailable, do not claim donor-robust state-level inference; report `INSUFFICIENT_DONOR_METADATA`.

### Acceptance checks

- Validate gene-ID mapping, genome build, ancestry and alleles, duplicate samples, sample and donor nesting, counts or normalized assay origin, per-state donor representation and missingness.
- Never treat cells from one donor as independent biological replicates.
- Flag state-specific donor dominance; e.g. a reference with 391 of 400 sampled cells from one donor has an effective source-donor count near 1, not 400 independent observations.
- State-specific data eligibility thresholds are fixed before genetic association results. Initial suggested **development gate:** at least five distinct donors per state **and** non-degenerate paired support. These are provisional; calibration and power evaluation may require higher thresholds.
- For case/control atlases, build labels using consistent covariate adjustment and compare matched biological contexts. A donor-based contrast is not equivalent to a cell-level Wilcoxon test.

## 3. Direction-aware annotation construction

For gene `g`, donor `d` and cell state `c`, define per-donor log-expression `x_gdc = log2(1 + CPM_gdc)` or use a valid donor-level normalized expression model (the choice is frozen). Only compare states that are supported **within the same donor and relevant strata**; never mix disease groups/phase unadjusted.

Define a paired donor contrast:

`Delta_gc = mean_d [ x_gdc - mean_{c' in matched comparators(d,c)} x_gdc' ]`.

Within a frozen, shared gene universe, construct the following candidate annotations:

- `UP_gc = max(Delta_gc, 0)` — gene enriched in the named state.
- `DOWN_gc = max(-Delta_gc, 0)` — gene depleted in the named state; **not** evidence that the GWAS acts in that state.
- `ABS_gc = abs(Delta_gc)` — a comparator mixing both directions.
- `SPEC_gc` — continuous cell-state specificity derived from matched donor-mean expression and normalized across reference states.

Provide both thresholded and predeclared fixed-size families (e.g. top 400 eligible genes) with **matched common universe** and diagnostics for set size, gene length, gene density, mappability, expression, and constraint. Ineligible/undersized sets receive `HOLD`, not zero enrichment.

**Interpretation rule:** UP/DOWN are directional *between-cell-type expression comparisons*. They do NOT encode direction of GWAS effect, eQTL effect, or disease-associated differential expression. A DOWN signal could indicate that risk genes act in a different cell state; its meaning needs direct cross-state checks.

## 4. Genetics engines and conditional tests

Run a frozen primary genetics engine with validated gene mapping, LD treatment and genomic baselines. For development, use MAGMA gene-property as primary and S-LDSC as secondary sensitivity. Reuse the official tools rather than approximate their p-values using independent-gene regressions.

Model standardized UP and DOWN annotations **jointly**, with their covariance and predefined generic-gene covariates; estimate `Delta` only if design rank, collinearity and available information meet quality gates. `SE(UP-DOWN)` must include the full covariance term. For related cells, condition focus annotation on competitors and report pairwise contrasts with multiplicity control. A nonsignificant difference does NOT establish equivalence or specificity.

Do not combine p-values from dependent annotations or overlapping GWAS as if they were independent studies. Standardize interpretations across MAGMA and S-LDSC; a first-ranked p-value is not sufficient evidence of the best cell.

## 5. Uncertainty axes and resolution policy

Perturb **only legitimate independent units**:

1. Donor leave-one-out / stratified donor bootstrap: rerun marker discovery and genetics, respecting donor/phase/batch structure.
2. LD-block jackknife or valid GWAS-aware resampling: measure genetic influence and covariance, respecting LD and sample overlap.
3. Annotation-rule sensitivity: UP, DOWN, ABS, SPEC, thresholds, gene universe, joint competitors.
4. Atlas sensitivity: independent atlas where available, with transparent overlap between reference studies.

Record rank stability, direction stability, effect and CI, group-vs-state consistency and sign/competitor overlap. **Bootstrap selection frequency is NOT a posterior probability.** Never turn these descriptive frequencies directly into '90% confidence the cell is causal.'

### Proposed hierarchical outputs

- `STATE_SUPPORTED`: prespecified multiplicity-adjusted evidence, no unresolved contrast with close competitors, passes donor and genetics eligibility gates, and a calibrated operating point on independent held-out benchmarks.
- `GROUP_SUPPORTED`: at least one lineage/group-level finding meets prespecified evidence and risk criteria, but specific states remain indistinguishable.
- `ABSTAIN`: inadequate evidence, failing stability/resolution, missing atlas states, or an uncovered biological setting.
- Administrative flags: `NOT_CALIBRATED`, `INSUFFICIENT_DONORS`, `GWAS_QC_HOLD`, `GENE_UNIVERSE_MISMATCH`, `REFERENCE_COVERAGE_HOLD`.

**Critical:** The actual risk thresholds and dominance statistics are **not set in v0.1**. Tune on training/validation data and lock before external test. Do NOT ship an apparently valid state/group assignment before calibration.

## 6. Calibration and reproducibility

### Internal development

- Separate `discovery` (Paper 9 data and earlier simulations), `design/tuning`, `calibration`, and `external test` by GWAS trait/family and atlas provenance. Shared genetic participants or embedded atlas subsets cannot constitute independent external validation.
- To report selective risk control, use truly held-out calibration/evaluation tasks, separate threshold selection from evaluation, and calculate binomial or other suitable confidence intervals for **incorrect assignments among answered cases**. Sparse calibration sets imply `NOT_CALIBRATED`.
- Provide risk–coverage and group-vs-exact-state curves; assess harm from abstention when missing cells are unrepresented by the reference.
- Calibrate statistical null error separately from biological **mislocalization** error: null FWER ≠ rate of choosing an incorrect state under a non-null GWAS.

### Benchmarks with known truth

(A) Null, (B) UP-driven truth, (C) DOWN-driven truth, (D) both-direction/shared signal, (E) cell identity/label ambiguity, (F) missing true cell type, (G) cell/donor imbalance, (H) ancestry/LD mismatch. Prefer LD-aware variant-level simulations or orthogonal regulatory truths that do not use the same RNA marker score to plant the source. Existing profile-rotation results remain **development evidence only**, not external validation.

### Baselines and primary metrics

Compare against (i) MAGMA/CELLECT specificity ranking, (ii) thresholded ABS and UP MAGMA/S-LDSC, (iii) scDRS, (iv) CATCH where runnable and transparent, and (v) a simple ensemble/abstention baseline. Report accuracy on **answered** tasks, state/group misassignment rates, risk–coverage AUC, null FWER/FDR, group-only yield, compute and failure rates. Match inputs and test same held-out tasks, not convenient subsets.

## 7. Versioning, scientific boundaries and stopping rules

- Paper 9 v001–v071 and prior submitted manuscript are **unchanged historical discovery results**. Do not retroactively change original verdicts or claim that development hypotheses were pre-registered before seeing prior outputs.
- P9-DARC starts as an **internal method-development project** with its own namespace and version control. No production-level or inferential claims from the toy prototype.
- Paper 7 retains PGR time-course and delay/attenuation identifiability. Paper 8 retains RNA/CpG markers, menstrual-fluid biology, and clinical diagnostics. Paper 9 uses RNA as a reference to interpret genetic enrichment, not as a new clinical biomarker panel.
- Stop and publish a negative methodological finding if held-out selective risk does not improve over baselines. Do not add traits/atlases only until a desired result appears.
- Functional ATAC/QTL follow-up is **optional biological corroboration**, not required for the method to run; orthogonal data must be held out from method tuning if used as external evidence.

## 8. v0.1 implementation status

The bundled `p9darc` prototype **only** validates small donor-level expression tables, computes directional donor-matched profiles, generates a transparent **noninferential** rank diagnostic against supplied MAGMA gene z-scores, and runs a LODO demonstration. Its output is ALWAYS `NOT_CALIBRATED`; it does **not** run MAGMA, S-LDSC, CATCH, selective inference or causal localization. Those are research milestones, not completed features.

## 9. Core references

- CELLECT: https://github.com/perslab/CELLECT
- scDRS: https://www.nature.com/articles/s41588-022-01167-z
- CATCH benchmark (Li et al.): https://www.medrxiv.org/content/10.1101/2025.05.24.25328275v3
- Broad prior-method survey: https://pmc.ncbi.nlm.nih.gov/articles/PMC11655331/
