import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# 1. PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures/Fig 3'))
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# High-impact aesthetics
sc.set_figure_params(dpi=150, dpi_save=600, format='pdf', transparent=True)

def main():
    print("--- PIPELINE SCRIPT 05: PAGA TRAJECTORY (FIGURE 3) ---")
    
    # 1. Load Discovery Cohort
    file_path = os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad')
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Missing discovery dataset: {file_path}")
        
    print("1. Loading Annotated Data...")
    adata = sc.read_h5ad(file_path)
    
    # Isolate Epithelial/Tumor cells for trajectory
    tumor_mask = adata.obs['cell_type'].astype(str).str.contains('Epithelial|Tumor|Cancer', case=False, na=False)
    adata_tumor = adata[tumor_mask].copy()
    
    # Ensure neighborhood graph exists
    if 'neighbors' not in adata_tumor.uns:
        print("Computing neighborhood graph...")
        sc.pp.neighbors(adata_tumor, n_neighbors=15, n_pcs=40, random_state=42)
        
    # 2. Run PAGA
    print("2. Computing PAGA topology (Root-free evolutionary lineage)...")
    
    # Define grouping column (ensure it matches the dataset annotations)
    group_col = 'Resistance_Status' if 'Resistance_Status' in adata_tumor.obs.columns else 'Sample_Origin'
    
    # Ensure categories exist
    if adata_tumor.obs[group_col].dtype.name != 'category':
        adata_tumor.obs[group_col] = adata_tumor.obs[group_col].astype('category')
        
    sc.tl.paga(adata_tumor, groups=group_col)
    
    # 3. Export PAGA Connectivity Matrix (Supplementary Data 3)
    print("3. Exporting PAGA Connectivity Matrix...")
    connectivities = adata_tumor.uns['paga']['connectivities'].todense()
    groups = adata_tumor.obs[group_col].cat.categories
    
    paga_df = pd.DataFrame(connectivities, index=groups, columns=groups)
    paga_csv = os.path.join(DATA_DIR, 'Supplementary_Data_3_PAGA_Connectivity.csv')
    paga_df.to_csv(paga_csv)
    print(f"   -> Saved: {paga_csv}")
    
    # 4. Plot Figure 3
    print("4. Generating Figure 3 (PAGA Graph)...")
    fig3_path = os.path.join(FIG_DIR, 'Figure_3_PAGA_Network.pdf')
    
    # Draw PAGA plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.paga(
        adata_tumor, 
        color=group_col, 
        edge_width_scale=1.5, 
        node_size_scale=2.0, 
        fontsize=10, 
        show=False,
        ax=ax,
        random_state=42
    )
    
    plt.savefig(fig3_path, bbox_inches='tight')
    plt.close()
    print(f"   -> Saved: {fig3_path}")
    
    print("--- SCRIPT 05 COMPLETE ---")

if __name__ == "__main__":
    main()
