# Osimertinib Resistance in NSCLC: Multi-Modal Transcriptomic Pipeline

This repository contains the computational pipeline used to investigate the transcriptomic mechanisms of Osimertinib resistance in Non-Small Cell Lung Cancer (NSCLC). 

## Overview
By integrating single-cell RNA sequencing (scRNA-seq), spatial transcriptomics, and machine learning, this framework identifies predictive genomic biomarkers and maps resistance evolution across diverse structural states. The pipeline highlights PCSK1N and ABCA3 as robust predictors of resistance, while structurally mapping a localized epigenetic seesaw between EZH2 and the HOXC genomic locus. It also demonstrates via in silico pseudo-bulking how traditional bulk sequencing methodologies mathematically mask these spatial correlations through Simpson’s Paradox.

---

## Repository Contents & File I/O Mapping
Run the scripts in the sequential order provided below. Below is the precise input and expected output mapping for each analytical step:

### `01_scRNAseq_preprocess_and_annotate.py`
Preprocessing & Quality Control: Ingests raw single-cell matrices, applies strict biological filtering (mitochondrial/doublet thresholds), normalizes data, and calculates high-dimensional embeddings (PCA/UMAP). Performs unsupervised Leiden clustering and annotates immune/tumor cell populations.
* **Input File:** `../data/GSE131907_RAW_merged.h5ad`
* **Output Files:** 
  * `../data/marker_genes_complete.csv`
  * `../data/GSE131907_PHASE1_ANNOTATED.h5ad`

### `02_predictive_modeling_and_biomarkers.py`
Biomarker Discovery: Isolates the epithelial tumor fraction to perform Differential Gene Expression (DGE) between Primary (sensitive) and Metastatic (resistant) states. Trains a Random Forest classifier to identify top predictive features (e.g., PCSK1N) and computes initial evolutionary pseudotime.
* **Input File:** `../data/GSE131907_PHASE1_ANNOTATED.h5ad`
* **Output Files:** 
  * `../data/Supplementary_Data_1_DGE_Cleaned.csv`
  * `../data/Supplementary_Data_2_Pathways.csv`
  * `../data/GSE131907_PHASE2_TUMORS.h5ad`

### `03_networks_and_clinical_validation.py`
Network Inference & Clinical Survival: Constructs the Protein-Protein Interaction (PPI) and Gene Regulatory Networks (GRN) for the identified biomarkers. Validates the clinical relevance of the signature using TCGA-LUAD patient cohorts via Kaplan-Meier survival curves and Multivariate Cox Proportional Hazards regression.
* **Input File:** `../data/TCGA_LUAD_Signature_Scores.csv`
* **Output Files:** 
  * `../data/Supplementary_Data_2_PPI_Edges.csv`
  * `../data/Network_Centrality_Analysis.csv`

### `04_spatial_transcriptomics_pipeline.py`
Spatial Niche Mapping: Processes 10x Visium spatial transcriptomics (GSE267960) to map the physical architecture of the resistance module. Tests for cell-cycle independence and maps the "Escape from Quiescence" spatial phenomenon in vivo.
* **Input File:** `../data/validation_cohorts/GSE267960/GSM8282531_S20-12254-L9/S20-12254-L9/filtered_feature_bc_matrix.h5`
* **Output Files:** 
  * `../data/GSE267960_Processed_Spatial.h5ad`
  * Console metrics displaying spatial activity and Pearson correlations (R-values) for cell cycle independence

### `05_functional_validation_and_epigenetics.py`
Epigenetics & Functional Genomics: Analyzes public genome-wide CRISPR-Cas9 synthetic lethality screens. Analyzes 3D organoid expression data to resolve in vitro 2D culturing artifacts. Simulates pseudo-bulk spatial masking to demonstrate Simpson's Paradox in standard RNA-seq workflows.
* **Input Files:** 
  * `../data/Pfeifer_et_al_Supplemental_Data/MAGECK_files/PC9_Osimertinib_CRISPRn.txt`
  * `../data/validation_cohorts/GSE267960/GSM8282531_S20-12254-L9/S20-12254-L9/filtered_feature_bc_matrix.h5`
* **Output Files:** 
  * `../data/CRISPR_Functional_Validation.csv`
  * `../data/Pseudo_Bulk_Masking_Simulation.csv`

### `06_robustness_and_confounder_ablation.py`
Advanced Rigor: Addresses reviewer critiques by executing a Stratified 5-Fold Cross-Validation on an ablated dataset (removing brain metastasis confounders). Generates unsupervised PAGA (Partition-based Graph Abstraction) evolutionary trajectories to prove root-independent lineage flow.
* **Input File:** `../data/GSE131907_FINAL_ANNOTATED.h5ad`
* **Output Files:** 
  * `../data/Cross_Validation_Metrics.csv`
  * `../data/Ablated_Feature_Importances.csv`
  * `../data/PAGA_Evolutionary_Connectivity.csv`

---

## Installation & Requirements
It is recommended to run this pipeline in an isolated Conda or Virtual Environment.

```bash
# Create a new conda environment
conda create -n osi-resistance python=3.9
conda activate osi-resistance

# Install required packages
pip install -r requirements.txt
```

**Key Dependencies:**
* `scanpy >= 1.9.0`
* `pandas >= 1.4.0`
* `numpy >= 1.21.0`
* `scikit-learn >= 1.0.0`
* `lifelines >= 0.27.0`
* `networkx >= 2.8.0`
* `gseapy >= 1.0.0`
* `scvelo >= 0.2.5`

---

## Data Availability
Due to GitHub file size limits, raw `.h5ad`, `.mtx`, and `.fastq` files are not included in this repository. Raw data can be downloaded directly from the Gene Expression Omnibus (GEO) and TCGA:

* **Main scRNA-seq Discovery Cohort:** GSE131907
* **3D Organoid Validation Cohort:** GSE255958
* **Spatial Transcriptomics Validation:** GSE267960
* **TCGA-LUAD:** Pan-Cancer Atlas 2018

Place all downloaded raw matrices into the `data/` directory before initiating Script 01. Generated for the manuscript resubmission regarding the landscape of Osimertinib resistance.