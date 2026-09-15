import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures/Fig 4'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 06: CELL-CELL COMMUNICATION (SQUIDPY) ---")
    
    # 1. Load Discovery Cohort
    adata_path = os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad')
    adata = sc.read_h5ad(adata_path)
    
    # Clean cluster annotations for permutations
    valid_cells_mask = ~adata.obs['cell_type'].astype(str).str.startswith('Cluster')
    adata_clean = adata[valid_cells_mask].copy()
    
    # Map tumor cells to Primary vs Metastatic
    def group_status(row):
        if row['cell_type'] == 'Epithelial cells':
            return 'Primary Tumor' if row['Sample_Origin'] in ['tLung', 'tL/B'] else 'Metastatic Tumor'
        return row['cell_type']
    
    adata_clean.obs['Detailed_Type'] = adata_clean.obs.apply(group_status, axis=1).astype('category')

    # 2. Run Permutations
    print("Running 1,000 Monte Carlo Permutations...")
    res = sq.gr.ligrec(
        adata_clean, cluster_key='Detailed_Type', n_perms=1000, 
        threshold=0.01, copy=True, seed=42
    )

    # 3. Export Statistical Matrices
    pvals_csv = os.path.join(DATA_DIR, 'Supplementary_Table_S4_Permutation_Pvalues.csv')
    means_csv = os.path.join(DATA_DIR, 'Supplementary_Table_S5_Permutation_Means.csv')
    res['pvalues'].to_csv(pvals_csv)
    res['means'].to_csv(means_csv)
    
    # 4. Generate Figure 4 (Dotplot of Key Axes)
    print("Generating Figure 4...")
    target_groups = ['Macrophages', 'T-Cells']
    source_groups = ['Primary Tumor', 'Metastatic Tumor']
    
    fig, ax = plt.subplots(figsize=(6, 4))
    sq.pl.ligrec(
        res, cluster_key='Detailed_Type', source_groups=source_groups, target_groups=target_groups,
        pvalue_threshold=0.05, alpha=0.05, swap_axes=True, ax=ax,
        save=os.path.join(FIG_DIR, "Figure_4_Immune_Evasion_Final.pdf")
    )
    print("--- SCRIPT 06 COMPLETE ---")

if __name__ == "__main__":
    main()
