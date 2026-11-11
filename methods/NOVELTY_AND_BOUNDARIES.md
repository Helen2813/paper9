# Literature overlap and separation from other papers (2026-11-11)

| Prior method | What it already does | What P9-DARC still needs to show beyond it |
| --- | --- | --- |
| CELLECT / CELLEX | Cell-type expression-specificity and MAGMA/S-LDSC genetics enrichment; multiple expression specificity metrics | Donor-induced uncertainty and *calibrated, hierarchical* refusal rather than another specificity rank |
| scDRS | Gene-based disease relevance at single-cell resolution with controls and cell-type aggregates | Directional-reference sensitivity and a model-independent state-vs-group resolution rule |
| CATCH / Li et al. (2026 preprint) | Extensive benchmarking and robust Cauchy combination of complementary prioritizers, reference sensitivity | Demonstrated out-of-sample selective risk reduction, not just a new score or modest rank shifts |
| Paper 9 prior v001–v071 | ABS/UP/DOWN/SPEC comparison, donor and locus audits, endometriosis biology | **Development data only**; cannot re-use as an unbiased external benchmark |

References: https://github.com/perslab/CELLECT ; https://www.nature.com/articles/s41588-022-01167-z ; https://www.medrxiv.org/content/10.1101/2025.05.24.25328275v3

## Portfolio boundaries

**Paper 7:** PGR time-series transfer/delay/attenuation/abstention detectors. No temporal PGR biology in P9-DARC. Conceptual similarity of abstention alone does not imply duplicated analysis.

**Paper 8:** biomarker transportability, RNA quintet TIMP2/AKR1C2/DMBT1/KCNK5/IL1A, CpG, menstrual effluent and clinical diagnostics. No new diagnostic panel claims in Paper 9.

**Paper 9:** inherited GWAS association, cell reference annotation, genomic enrichment, cell-attribution uncertainty and external method benchmarks. Sharing atlases is permissible; reuse of the same cases/genes cannot be called independent validation.

## BioRxiv and journal positioning

An implemented tool is not automatically original research. A future Method/Software manuscript must show validated novelty against prior methods. The submitted Paper 9 and its bioRxiv screening decision are historical records, not erased by this new project. Do not promise acceptance solely from rewriting or releasing a CLI.
