import os
import sys
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 1'))
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

sc.set_figure_params(dpi=150, dpi_save=600, format='pdf', transparent=True, color_map='viridis')

def main():
    print("--- PIPELINE SCRIPT 01: GSE131907 PREPROCESSING ---")
    raw_path = os.path.join(DATA_DIR, 'GSE131907_RAW_merged.h5ad')
    if not os.path.exists(raw_path):
        print(f"[ERROR] Missing raw data: {raw_path}")
        sys.exit(1)

    adata = sc.read_h5ad(raw_path)
    if not sp.issparse(adata.X): adata.X = sp.csr_matrix(adata.X)
    
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[(adata.obs.n_genes_by_counts < 5000) & (adata.obs.pct_counts_mt < 10), :].copy()
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=40, random_state=42)
    sc.tl.umap(adata_hvg, random_state=42)
    sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)
    
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    
    cluster_annotations = {'1': 'T-Cells', '2': 'Macrophages', '3': 'NK / Cytotoxic Cells', '4': 'B-Cells'}
    for cluster in ['9', '11', '12', '13', '20', '21', '24', '25', '26']: cluster_annotations[cluster] = 'Epithelial cells'
    adata.obs['cell_type'] = adata.obs['leiden'].map(lambda x: cluster_annotations.get(x, f"Cluster {x}")).astype('category')
    
    origin_mapping = {'nLung': 'Normal', 'nLN': 'Normal', 'tLung': 'Primary', 'tL/B': 'Primary', 'mLN': 'Metastatic', 'PE': 'Metastatic', 'mBrain': 'Metastatic'}
    adata.obs['Status'] = adata.obs['Sample_Origin'].map(origin_mapping)
    
    fig1a = sc.pl.umap(adata, color='cell_type', palette=sns.color_palette("colorblind").as_hex(), show=False, return_fig=True)
    fig1a.savefig(os.path.join(FIG_DIR, 'Figure_1A_CellTypes.pdf'), bbox_inches='tight')
    plt.close(fig1a)
    
    fig1b = sc.pl.umap(adata, color='Status', palette=sns.color_palette("colorblind").as_hex(), show=False, return_fig=True)
    fig1b.savefig(os.path.join(FIG_DIR, 'Figure_1B_ClinicalOrigin.pdf'), bbox_inches='tight')
    plt.close(fig1b)

    adata.write(os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad'))
    print("--- SCRIPT 01 COMPLETE ---")

if __name__ == "__main__":
    main()
