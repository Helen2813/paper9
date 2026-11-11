# P9-DARC — Benchmark and Validation Plan v0.1

**Status:** Prespecified *proposal* for a new method. Full calibration design and task cohorts must be frozen before first external benchmark. The historical endometriosis results are not held out.

## Experimental table

| Arm | Test | Required independent truth / control | Fail criterion |
| --- | --- | --- | --- |
| A: Null | Null GWAS and LD-preserving resampling | No planted genetics | Empirical family-wise null error above prespecified tolerance, or malformed null |
| B: Direction | Genuine UP, DOWN, mixed latent variant mechanisms | Independent causal variants / known regulatory mapping | Bad selective risk or inferior accuracy at matched coverage |
| C: Donors | 1 large donor vs many balanced donors, matched donor counts | Held-out donor identity | Spurious confidence rises after donor imbalance |
| D: Missing state | True state removed from reference; parent group may remain | Independent held-out state identity | Wrong sibling assigned without adequate abstention |
| E: Related states | Similar fibroblast/perivascular/immune subclasses | Frozen hierarchy from independent annotations | Unjustified fine-state calls; wrong group |
| F: Atlases | Different platform / institution / assay references | Independently generated atlas | Serious unreported transport failure |
| G: External traits | Endocrine, immune, vascular, neurological with valid atlas coverage | Predeclared broad biology, not presumed single causal cell | No improvement over strong baseline on held-out tasks |
| H: Genetic | Different LD architectures and GWAS powers | Variant-level truth | Inflated claims at low power or ancestry mismatch |

## Evaluation requirements

1. **Freeze** candidate architectures, hierarchy definitions, genetic engine, null construction, thresholds, and metric definitions on *development* tasks.
2. Partition at **trait-family/GWAS provenance level** and **atlas-study level**; endometriosis (Koller / GCST / FinnGen) is a discovery example, NOT three independent external tasks.
3. Use at least two non-endometriosis biological systems with appropriate independent single-cell atlases in test. Exact number depends on feasible data; if too few independent truths, report descriptive results without calibrated generalization.
4. Track `STATE_SUPPORTED`, `GROUP_SUPPORTED`, `ABSTAIN`, `NOT_CALIBRATED` separately; score exact-state and parent-group errors on the answered subset and across all cases.
5. Plot selective risk (error among predictions) vs coverage; compare to naive max-p/lowest-p, UP-only, CELLECT, scDRS and CATCH when execution is feasible. Method cannot claim superiority if baseline comparison is absent.
6. Provide CIs for *paired* contrasts via correct grouping by independent tasks. Gene-wise bootstrap is not a valid substitute for LD structure or independent traits.
7. Report positive and negative results, all failed runs, gene / SNP / donor attrition, ancestry/build metadata, versions and hashes.
8. If no state-level assignments pass calibration but group assignments do, that is an acceptable method result. If no improved risk–coverage is demonstrated, report the negative benchmark and do not call it a universal solution.

## Claim gate matrix

| Claim | Evidence required | Current status |
| --- | --- | --- |
| Direction changes rankings | Existing Paper 9 analyses | **Observed in development data** |
| Donor balance matters | Paper 9 v054 | **Observed in development data** |
| Abstention is calibrated | Held-out truth and independent selection/evaluation | **NOT TESTED** |
| Better than CELLECT/CATCH | Matched input benchmark | **NOT TESTED** |
| Generalizes outside endometriosis | Several independent trait/atlas pairs | **NOT TESTED** |
| Causal cell identified | Valid causal functional genetics | **NOT ESTABLISHED; outside primary target** |

## Stage gates and scope

- `S0` (current): method specification and an audit-only prototype + unit tests.
- `S1`: exact input loaders for v054/v055 (read-only) and official MAGMA/S-LDSC wrappers with reference provenance; no selection rule tuning.
- `S2`: valid donor/LD resampling and known-truth variant-level simulation, including missing states.
- `S3`: internal calibration and separate external evaluation, matched comparator runs.
- `S4`: finalized R package / Python CLI only if S3 shows reproducible benefit; otherwise publish an evidence-limited methodological audit.

**Note:** A positive ATAC result in GSE183771 can corroborate the biology of the endometriosis case study but does not substitute for calibrated external benchmarking of P9-DARC.
