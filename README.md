This repository contains the core computational pipelines and analysis scripts for Ancestry-Associated Molecular and Cellular Immunosuppresive Features at Single Cell
Resolution in Treatment-Naive Non-Small Cell Lung Cancer.

```text
NSCLC-scRNA-seq/
├── .gitignore
├── README.md
├── session_info_main_seurat.txt         # Package versions for Seurat/CellChat workflow
├── session_info_trajectory_monocle2.txt # Package versions for Monocle2 trajectory workflow
├── session_info_supp_glmmtmb.txt        # Package versions for general linear mixed model workflow
└── scripts/
    ├── seurat_main/
    │   ├── fig_01_global_characteristics.R
    │   ├── fig_02_epithelial_compartment.R
    │   ├── fig_03_stromal_compartment.R
    │   ├── fig_04_myeloid_compartment.R
    │   ├── fig_05_t_compartment.R
    │   ├── fig_06_b_compartment.R
    │   ├── fig_07_cell_cell_interactions.R
    │   ├── supp_01_global_characteristics.R
    │   ├── supp_02_epithelial_compartment.R
    │   ├── supp_03_stromal_compartment.R
    │   ├── supp_04_myeloid_compartment.R
    │   ├── supp_05_t_compartment.R
    │   ├── supp_06_b_compartment.R
    │   └── supp_07_cell_cell_interactions.R
    │
    ├── glmm_patient_level/
    │   ├── Supp_module_scores_GLMM.R
    │   ├── Supp_SPP1_GLMM_sensitivity.R
    │   ├── Supp_TAM_SPP1_subclustering.R
    │
    └── monocle2_trajectory/
        ├── monocle2_pseudotime_myeloid.R
        ├── monocle2_pseudotime_cd8t.R
        ├── monocle2_pseudotime_cd4t.R
        └── monocle2_pseudotime_epithelial.R
```

The raw and processed single-cell transcriptomic datasets for this study are deposited at [GEO:]. The scripts in this repository are organized by figure panels and are intended for analytical review, parameter checking, and reproducibility auditing. Users wishing to re-run specific pipelines locally will need to download the source data objects and update the setwd() or file paths at the top of each script to match their local storage environment.


| Script Name | Target Manuscript Panels | Description / Analysis Focus |
| :--- | :--- | :--- |
| `fig_01_global_characteristics.R` | Figure 1 | Cellular heterogeneity in the TME of Black and White NSCLC patients |
| `fig_02_epithelial_compartment.R` | Figure 2 | A Hyper-Metabolic, Basal-like Malignant Epithelial Program |
| `fig_03_stromal_compartment.R` | Figure 3 | Stromal cell profiling (Fibroblasts, Endothelial) |
| `fig_04_myeloid_compartment.R` | Figure 4 | Myeloid cell profiling (Macrophages, Monocytes, DCs) |
| `fig_05_t_compartment.R` | Figure 5 | T cell subset heterogeneity and functionality |
| `fig_06_b_compartment.R` | Figure 6 | B cell and Plasma cell subset |
| `fig_07_cell_cell_interactions.R` | Figure 7 | Global network communication |
| `supp_01_global_characteristics.R` | Supplemental Figure 1 | Supplemental Global Landscape |
| `supp_02_epithelial_compartment.R` | Supplemental Figure 2 | A Hyper-Metabolic, Basal-like Malignant Epithelial Program |
| `supp_03_stromal_compartment.R` | Supplemental Figure 3 | Stromal cell profiling (Fibroblasts, Endothelial) |
| `supp_04_myeloid_compartment.R` | Supplemental Figure 4 | Myeloid cell profiling (Macrophages, Monocytes, DCs) |
| `supp_05_t_compartment.R` | Supplemental Figure 5 | T cell subset heterogeneity and functionality |
| `supp_06_b_compartment.R` | Supplemental Figure 6 | B cell and Plasma cell subset |
| `supp_07_cell_cell_interactions.R` | Supplemental Figure 7 | Comprehensive CellChat panel outputs |
| `Supp_module_scores_GLMM.R` | Supplemental Tables | GLMM analysis for modular scores |
| `Supp_SPP1_GLMM_sensitivity.R` | Supplemental Tables | GLMM analysis for SPP1 expression |
| `Supp_TAM_SPP1_subclustering.R` | Supplemental Tables | SGLMM analysis for SPP1 expression in TAM subclusters |


Contact: For questions regarding the code or data access, please contact [Pengbo Zhang] at [pengbo.zhang@advocatehealth.org].

The repository was updated during peer review to use two-sided Wilcoxon rank-sum tests with Benjamini–Hochberg FDR correction and to include the supplementary GLMM analyses for SPP1 expression and modular scores.

Citation: [...]
