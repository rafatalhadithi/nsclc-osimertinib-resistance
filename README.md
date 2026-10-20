# Osimertinib Resistance in NSCLC: Multi-Modal Transcriptomic Pipeline

This repository contains the complete, reproducible computational pipeline used to investigate the transcriptomic mechanisms of Osimertinib resistance and metastatic lineage plasticity in Non-Small Cell Lung Cancer (NSCLC). 

## Overview
By integrating single-cell RNA sequencing (scRNA-seq), spatial transcriptomics, and machine learning, this framework maps resistance evolution across diverse structural states. A Nested Patient-Level Cross-Validated Random Forest classifier (GroupKFold) highlights **SCGB3A2** and **MUC1** as robust drivers of acquired resistance. The pipeline utilizes `decoupler` for rigorous pseudo-bulked GSEA, and `squidpy` permutation testing to identify the dynamic acquisition of the **LGALS9-HAVCR2 (TIM-3)** immune exhaustion axis during Progressive Disease. Finally, targets are validated across TCGA clinical cohorts, revealing that tumors actively hijack indolent baseline epithelial states to survive acute targeted therapy.

---

## Data Availability
Due to GitHub file size limits, raw `.fastq`, `.mtx`, and massive `.h5ad` single-cell matrices are not included in this repository. 

**For Peer Reviewers (Pre-Processed Data):**
To facilitate immediate reproducibility without requiring extensive computational preprocessing, the fully processed and annotated single-cell AnnData objects used in this study are securely hosted on Zenodo.
1. Download `GSE131907_FINAL_ANNOTATED.h5ad` (Discovery) and `maynard_annotated.h5ad` (Validation) from Zenodo: **[https://doi.org/10.5281/zenodo.20999292](https://doi.org/10.5281/zenodo.20999292)**
2.  Place the downloaded `.h5ad` files into the `data/` directory of this cloned repository before running the scripts.

**Raw Data Sources:**
Raw matrices can be accessed via the Gene Expression Omnibus (GEO) and Sequence Read Archive (SRA):
* **Discovery Cohort:** GSE131907
* **Longitudinal Validation Cohort:** PRJNA591860 (Maynard et al., 2020)
* **Spatial Transcriptomics Validation:** GSE267960
* **TCGA-LUAD:** Pan-Cancer Atlas 2018

---

## File Name Mapping Guide
All generated data tables have been standardized to map directly to the manuscript's Supplementary Data.

| Manuscript File Name | Description |
| :--- | :--- |
| `Supplementary_Data_1.csv` | Pseudo-bulked differential gene expression (DGE) results. |
| `Supplementary_Data_2.csv` | MSigDB Hallmark GSEA enrichment results (Discovery Cohort). |
| `Supplementary_Data_3.csv` | PAGA evolutionary trajectory and connectivity matrix. |
| `Supplementary_Data_4.csv` | GroupKFold CV performance metrics (Accuracy, F1, Precision, Recall). |
| `Supplementary_Data_5.csv` | Top biomarker genes (SCGB3A2/MUC1) and RF importance weights. |
| `Supplementary_Data_6.csv` | Spot-level Spearman correlations between epigenetic regulators (EZH2-HOXC). |
| `Supplementary_Data_7.csv` | Evaluation of spatial masking via Simpson's Paradox (Simulated Bulk). |
| `Supplementary_Data_8.csv` | EZH2 cross-cohort architectural variance metrics. |
| `Supplementary_Data_9.csv` | Maynard single-cell validation data (Expression metrics for Violin Plots). |
| `Supplementary_Data_10.csv` | Cross-cohort GSEA validation summary (Maynard PD vs TN). |
| `Supplementary_Data_11.csv` | Maynard patient-level metadata manifest and retained cell counts. |

---

## Pipeline Execution Details
All scripts utilize dynamic relative paths for cross-platform compatibility. Execute the scripts from within the `notebooks/` directory in the sequential order provided below:

### `01_preprocess_GSE131907.py`
**Preprocessing & Quality Control:** Ingests raw matrices, applies strict biological filtering, normalizes data, and calculates high-dimensional embeddings (PCA/UMAP/Leiden). Annotates populations as Primary or Metastatic.
* **Inputs:** `GSE131907_RAW_merged.h5ad`
* **Outputs:** `GSE131907_FINAL_ANNOTATED.h5ad` [Zenodo], Figure 1, Supp Fig S1, Supp Fig S2

### `02_preprocess_Maynard.py`
**Validation Cohort Preparation:** Standardizes the external PRJNA591860 longitudinal dataset. Extracts clinical TN, RD, and PD stages, deduplicates metadata, and filters for high-quality tumor and immune cells.
* **Inputs:** `S01_metacells.csv`, `S01_datafinal.csv`
* **Outputs:** `maynard_annotated.h5ad` [Zenodo], `Supplementary_Data_11.csv`

### `03_ML_Patient_CV.py`
**Confounder-Free Biomarker Discovery:** Isolates the epithelial tumor fraction and utilizes a Nested Patient-Level GroupKFold CV Random Forest classifier to prevent data leakage. Identifies top predictive features.
* **Inputs:** `GSE131907_FINAL_ANNOTATED.h5ad`, `maynard_annotated.h5ad`
* **Outputs:** `Supplementary_Data_4.csv`, `Supplementary_Data_5.csv`, Figure 4

### `04_PseudoBulk_GSEA.py`
**Pathway Enrichment:** Uses `decoupler` to aggregate cells by patient to eliminate pseudo-replication inflation. Runs continuous Log2FC-ranked GSEA against MSigDB Hallmark signatures for both cohorts.
* **Inputs:** `GSE131907_FINAL_ANNOTATED.h5ad`, `maynard_annotated.h5ad`
* **Outputs:** `Supplementary_Data_1.csv`, `Supplementary_Data_2.csv`, `Supplementary_Data_10.csv`, Figure 2

### `05_PAGA_Trajectory.py`
**Evolutionary Lineage:** Generates unsupervised PAGA (Partition-based Graph Abstraction) evolutionary trajectories utilizing a Force-Directed layout to map structural continuums.
* **Inputs:** `GSE131907_FINAL_ANNOTATED.h5ad`
* **Outputs:** `Supplementary_Data_3.csv`, Figure 3

### `06_Squidpy_Maynard_Evasion.py`
**Cell-Cell Communication:** Replaces descriptive co-expression with a 1,000-permutation statistical framework via `squidpy` on the high-depth Smart-seq2 cohort to statistically validate CD274-PDCD1 and LGALS9-HAVCR2 (TIM-3) immune evasion interactions.
* **Inputs:** `maynard_annotated.h5ad`
* **Outputs:** Figure 4 (Note: Saved internally as Figure 4 Immune Evasion)

### `07_Maynard_Violins.py`
**Clinical Stage Validation:** Longitudinally tracks single-cell expression shifts of *SCGB3A2*, *MUC1*, *MUC5B*, and *EZH2* across TN, RD, and PD states.
* **Inputs:** `maynard_annotated.h5ad`
* **Outputs:** `Supplementary_Data_9.csv`, Figure 6

### `08_Spatial_Simpson_Paradox.py`
**Spatial Niche Mapping:** Simulates pseudo-bulk spatial masking of the *EZH2-HOXC* epigenetic seesaw to mathematically demonstrate Simpson's Paradox.
* **Inputs:** None (Statistical simulation executing extracted coordinates from *in vivo* GSE267960 Visium spatial data).
* **Outputs:** `Supplementary_Data_6.csv`, `Supplementary_Data_7.csv`, Supp Fig S3, Supp Fig S4

### `09_PPI_Network_Generator.py`
**Network Architecture:** Constructs the Protein-Protein Interaction (PPI) network for the identified ML biomarkers via the STRING API, generating a Cytoscape-ready topology.
* **Inputs:** `Supplementary_Data_5.csv`
* **Outputs:** Figure 8 (.graphml format for Cytoscape)

### `10_TCGA_Survival_Prep.py`
**Clinical Stratification Prep:** Parses TCGA-LUAD EGFR+ expression and clinical metadata, calculates Composite Z-Scores of the *MUC1/SCGB3A2* biomarker signature, and pre-processes quartile survival groups.
* **Inputs:** `Supplementary_Data_5.csv`, `data_mrna_seq_v2_rsem.txt`, `data_clinical_patient.txt`, `egfr_status.txt`
* **Outputs:** `Figure_7B_Survival_Data.csv`, `Figure_9_Survival_Data.csv`

### `11_Plot_Survival_R.R`
**Kaplan-Meier Survival Plotting:** R script utilizing `survminer` to calculate Cox Proportional Hazards and generate high-impact Progression-Free Survival (PFS) and Overall Survival (OS) curves with at-risk tables.
* **Inputs:** `Figure_7B_Survival_Data.csv`, `Figure_9_Survival_Data.csv`
* **Outputs:** Figure 7B, Figure 9

---

## Installation & Requirements
To ensure full reproducibility of the computational pipeline, the following package versions (or higher) are required. We strongly recommend running this pipeline within a dedicated Conda environment.

```bash
conda create -n osi-resistance python=3.13
conda activate osi-resistance
pip install scanpy pandas numpy scikit-learn decoupler-py squidpy gseapy networkx lifelines matplotlib seaborn


#  Key Dependencies

To ensure full reproducibility of the computational pipeline, the following package versions (or higher) are required. We strongly recommend running this pipeline within a dedicated virtual environment or Conda environment.

*   **scanpy** >= 1.9.0
*   **pandas** >= 1.4.0
*   **numpy** >= 1.21.0
*   **scikit-learn** >= 1.0.0
*   **decoupler** >= 1.4.0
*   **squidpy** >= 1.2.0
*   **gseapy** >= 1.0.0
*   **networkx** >= 2.8.0
*   **lifelines** >= 0.27.0

