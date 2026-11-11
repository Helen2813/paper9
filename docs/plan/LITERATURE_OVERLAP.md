# Paper 9 — literature and dataset overlap (v0.2.11, 2026-11-06)

## Datasets

| Dataset | Correct source | Used in Paper 9 | Inside HECA? |
|---|---|---|---|
| GSE213216 | Fonseca et al. 2023, Nat Genet 55:255-267 (Lawrenson lab) | G1 (lesion vs eutopic) | yes |
| GSE179640 | Tan et al. 2022, Nat Cell Biol | not used directly (its nomenclature labels GSE266265) | yes |
| GSE266265 | 60 donors (35 endometriosis, 25 controls), GEO 2026, ScaiVision study | G4 checks 1-2, G1 reference centroids | no |
| HECA | Mareckova et al. 2024, Nat Genet (Wang, Garcia-Alonso, Tan, Lai, Fonseca, Huang + new) | G4 check 3-5 | — |
| GCST90205183 | Rahmioglu et al. 2023, Nat Genet (excl. 23andMe; includes FinnGen) | all genetics | — |

Erratum: PLAN_v0.2 and the v015-v018 script headers call GSE213216 "Tan 2022"; this is wrong.

## Published work that overlaps Paper 9 questions

- Rahmioglu et al. 2023 (Nat Genet): TissueEnrich on genes near 49 index SNPs -> endometrium,
  uterus, smooth muscle. Only top loci; no polygenic cell-type test.
- Mareckova et al. 2024 (HECA, Nat Genet): fGWAS (+-500 kb TSS, mean expression per type) with
  GCST90205183 -> dStromal early/mid and macrophages. Paper 9 G4 tests its robustness.
- Guare et al. 2025 (GBMI; PMC12687807): 50,265 cases; scDRS on HECA comparing cells of donors
  with vs without endometriosis; strongest proliferative stroma. Case-control within eutopic
  endometrium, not heritability enrichment and not lesion vs eutopic.
- Ochoa et al. 2025 (Adv Sci; bioRxiv 2024.11.21.624726): scDRS on GSE213216 lesion cells only
  (200,798 cells; endometrioma + peritoneal), GWAS 23,492 cases; M2 macrophages p = 0.013,
  endothelial 0.044, dendritic 0.050 (uncorrected Monte Carlo p, no trait-specificity controls
  reported); IL1 assembloid and mouse experiments. No eutopic-cell scDRS result reported.
- Koller et al. 2026 (Nat Genet, doi 10.1038/s41588-026-02582-2): multi-ancestry, 105,869 cases,
  80 loci (37 novel); TWAS (S-PrediXcan / S-MultiXcan) in uterus, vagina, ovary, breast, whole
  blood, colon sigmoid / transverse (41 genes); colocalization; includes 23andMe.
  Published version (author manuscript, PMC 2026-05-22): S-LDSC on HECA marker-gene annotations
  (|log2FC| > 1.2, FDR < 0.01, cases and controls separately) -> "nominally significant
  enrichment for decidual stromal cells and SOX9 functionalis" (Suppl. Tables 9-10).
  Summary statistics: Zenodo doi 10.5281/zenodo.18983492; 23andMe part restricted.
  Implication for G4: with ~1.7x our cases, the decidual signal is again nominal only; it
  converges with our checks 1-3 (positive, high rank, not significant after correction) and
  with our check 4a (best HECA type SOX9_functionalis_I).

## Published cell-context claims to stress-test (Option A, PLAN_v0.3_DRAFT)

| Claim | Source | Method | GWAS |
|---|---|---|---|
| Decidual stroma (early/mid) + macrophages | HECA 2024 | fGWAS, +-500 kb TSS, mean expression | GCST90205183 |
| M2 macrophages (+ endothelial, DC) in lesions | Ochoa 2025 | scDRS, top 1,000 MAGMA genes | 23,492 cases |
| Proliferative stroma (case vs control donors) | Guare 2025 | scDRS on HECA | GBMI 50,265 cases |
| Decidual stroma, SOX9 functionalis (nominal only) | Koller 2026 | S-LDSC, HECA marker genes | 105,869 cases |

## What remains distinct in Paper 9

1. A pre-registered, multi-method, multi-atlas stress test of the published decidual-stroma
   enrichment (G4), including why methods disagree (check 4: not a few loci; check 5: annotation
   design).
2. A direct lesion-vs-eutopic contrast of cell programs against genetics (G1), with
   trait-specificity controls run for every contrast (an unrelated trait reproduced the only
   signal).
3. Explicit detection limits for every null result.
4. G2: directional concordance of genetically predicted expression with lesion-vs-eutopic versus
   patient-vs-control expression (not done by the TWAS papers).
5. Subtype architecture (G3) with FinnGen R13 subtypes.

## Added 2026-11-08
| Li et al. 2026 (medRxiv 10.1101/2025.05.24.25328275), "Benchmarking methods integrating GWAS and single-cell transcriptomic data for mapping trait-cell type associations" | 20 method configurations, simulated GWAS with known beta-cell source, FWER, misassignment to related types when the causal type is missing, CATCH | Must be cited; Paper 9 differs: directional marker selection (|log2FC| / up / down), a published endometriosis case, related endometrial states |
| Olafsson et al. 2026 (bioRxiv 10.64898/2026.02.25.708037), endometriosis lesions oligoclonal, derived from normal endometrium (WGS of 403 microdissections) | Cell of origin by somatic genomics | Context for discussion; Paper 9 does not address lesion origin |

## Added 2026-11-11 (chromatin test, v0.2.61)

- Garcia-Alonso, Vento-Tormo et al. 2026, bioRxiv 10.64898/2026.06.10.731198 ("An integrated multimodal pan-organ atlas of
  the female reproductive system across the lifespan"): CELLECT-LDSC on scRNA expression specificity -> endometriosis
  heritability in late-secretory stromal fibroblasts and smooth muscle; scATAC (43 donors, 718,670 peaks) used for
  peak-gene links and GWAS locus overlap, not for S-LDSC on cell-type peaks. Supports the stromal answer of Paper 9 by
  another expression-based method; must be cited. Peaks released only as a union list (GitHub ventolab).
- Vrljicak et al. 2023, Cell Rep ("Dynamic chromatin remodeling in cycling human endometrium at single-cell level";
  GSE183771): scATAC of 12 timed biopsies; no GWAS analysis. Source of the v072 lineage chromatin annotations.
- GSE246585 (EnSC decidualization ATAC time course): used only as an external check of stromal peaks (time points pooled).
