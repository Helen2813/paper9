# Paper 9 — Directional marker selection and genetic cell-state prioritization in endometriosis

**Research manuscript and reproducibility materials**  
**Authors:** Elena Krikun and Abdelrahman Alkhateeb  
**Status:** Manuscript in preparation / journal submission not represented by this repository. The associated bioRxiv submission was not posted; a repository is not a peer-reviewed publication.

## Overview

Genome-wide association studies (GWAS) can identify regions linked to endometriosis, but assigning those associations to individual cell states requires additional assumptions. A common strategy tests whether genes selected as single-cell markers are enriched for GWAS signal. **The direction used to select those markers can change which cell state ranks first.**

This study reconstructs published endometriosis cell-state prioritization, compares up-regulated, down-regulated, and direction-agnostic marker selection, evaluates calibration and trait controls, and asks whether any biological conclusion is robust to donors, reference atlases, and an independent molecular assay (single-nucleus ATAC-seq).

The central conclusion is **evidence for a broad stromal regulatory context for part of inherited endometriosis risk**, not proof of a uniquely causal stromal subtype or the cellular origin of endometriotic lesions.

## Main findings

1. **Marker direction changes the answer.** Direction-agnostic selection prioritized venous endothelium in a reconstructed analysis; up-regulated markers prioritized decidual stroma. The latter ranking did not survive within-family false-discovery correction in that comparison.
2. **No marker rule performed best across all checks.** A calibrated benchmark and a six-trait panel exposed sensitivity to how marker sets are defined and to related-cell-state ambiguity.
3. **Stromal interpretations are more stable than specific stromal labels.** Donor-level checks and a separate RNA atlas repeatedly implicated stromal programs, while exact fibroblast or perivascular assignments varied. Evidence for mesothelium remained inconsistent across cohorts.
4. **An orthogonal chromatin assay supported the compartment-level result.** In GSE183771, 36,073 nuclei from **12 biopsies** were retained. The initial single-annotation stromal S-LDSC test gave **z = 4.11**, one-sided **p = 2 × 10^-5**.
5. **Post-result joint sensitivity models retained a stromal signal.** With baselineLD v2.2 and lineage annotations entered together, KOLLER_EUR stromal enrichment was **z = 3.828**, one-sided **p = 6.471 × 10^-5**. The covariance-aware stroma-minus-epithelium contrast was **z = 3.677**, **p = 1.182 × 10^-4**.
6. **Biopsy-level sensitivity did not depend on one specimen.** In 12 leave-one-biopsy-out refits, stromal z statistics ranged from **3.439 to 3.890**, with positive nominal associations in every refit. These are correlated sensitivity tests, **not 12 independent replications**.

The joint model also found a stromal signal for **leiomyoma**, so this is not an endometriosis-exclusive genetic association. The endometriosis GWAS are partly overlapping and must not be counted as independent genetic replications.

## Files in the accompanying manuscript package

The Overleaf manuscript archive and this GitHub repository may have different contents. The following paths refer to the provided *manuscript package*; do not assume every item has already been committed to GitHub.

| Path | Purpose |
| --- | --- |
| `main.tex` | Main manuscript source (the current revised Overleaf version is self-contained) |
| `figures/` | Manuscript and supplementary figures, F1–F11 and S1–S7 |
| `sections/` | Section-level manuscript source files from the original multi-file package |
| `tables/` | Main/supplementary LaTeX tables and an editable Table 2 CSV |
| `code/p9_v074_scatac_hier.py` | Chromatin lineage-label/annotation reconstruction |
| `code/paper9_v076h.py` | Post-result joint chromatin S-LDSC sensitivity audit |
| `code/paper9_v077d.py` | Biopsy-exclusion annotation and S-LDSC audit |
| `data/v076h_joint_all_traits.csv` | Summary of joint-model estimates across six traits and two baselines |
| `data/v077d_biopsy_loo.csv` | Summary of the 12 biopsy-exclusion sensitivity estimates |
| `MANIFEST_SHA256.json` | Checksum manifest for the original manuscript export |

**Important:** `MANIFEST_SHA256.json` describes the original export. If you replace `main.tex` or add this README, do not claim that the old manifest verifies the modified project; regenerate it if checksums are needed.

## Reproducibility and provenance

Analyses use public GWAS summary statistics and public single-cell resources, including HECA, GSE266265, GSE213216, GSE179640, and GSE183771. See the manuscript's Methods and References for accession-level and method citations. Relevant methods include MAGMA and stratified LD score regression (S-LDSC).

The provided scripts document specific *analysis stages*, not an independently executable one-command end-to-end workflow. In particular:

- Some scripts assume a local directory at `~/Desktop/paper9`, cached input files, an existing project virtual environment, installed scientific software, and previously generated LD scores or annotations.
- `paper9_v077d.py` depends on earlier `v076h` code and frozen `v074` sources, which must be available in compatible locations; the included summary CSV is not a substitute for those inputs.
- Raw data, full intermediate outputs, and all earlier analysis scripts are **not included** in the manuscript archive. Their absence limits standalone reproduction from this snapshot.
- Data licenses, access terms, sample-level provenance and applicable redistribution restrictions must be checked before mirroring inputs.
- Initial hypothesis tests and post-result sensitivity analyses are distinguished in the manuscript; please preserve that separation when reusing results.

Selected audit outputs can be inspected directly in:

```text
data/v076h_joint_all_traits.csv
data/v077d_biopsy_loo.csv
```

**Do not infer clinical diagnostic utility, patient-level predictive performance, a unique causal cell state, or lesion cell of origin from these enrichment results.**

## Manuscript

The LaTeX article contains the full study, methods, figures and supplementary material. For the current single-file revision, upload the revised `main.tex` alongside its referenced figure files and compile using **XeLaTeX** (the `fontspec` package is incompatible with pdfLaTeX).

The figures are intended to appear near their first substantive discussion in the manuscript rather than as one block at the end. Supplementary figures and tables remain identified separately.

## Citation

Please cite the published article or publicly posted preprint **once one exists**. Until then, refer to this work as an unpublished manuscript and include the authors, title, and version/date accessed as appropriate. **No DOI is claimed for this repository or manuscript here.**

## Scope and limitations

This is a methods and biological-context study, not a validated clinical test. Results are sensitive to marker choice and reference atlas; stromal/perivascular labels overlap; some datasets are small or donor-imbalanced; the chromatin biopsies are not an endometriosis case–control cohort; and several GWAS have overlapping participants. The chromatin cell labels also use RNA marker information, so that analysis is independent of RNA **marker-set selection**, not completely marker-free.

## Contact and reuse

For questions about data provenance or reproducing an analysis, please open a GitHub issue in the repository (if issues are enabled). No open-source license is asserted by this README; reuse rights are governed by the repository's actual license, if any, and the upstream datasets' terms.
