import os
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 1. PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures'))

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(os.path.join(FIG_DIR, 'Fig 1'), exist_ok=True)
os.makedirs(os.path.join(FIG_DIR, 'Supplementary Figures', 'S1'), exist_ok=True)
os.makedirs(os.path.join(FIG_DIR, 'Supplementary Figures', 'S2'), exist_ok=True)

# Aesthetics
sc.set_figure_params(dpi=150, dpi_save=600, format='pdf', transparent=True, color_map='viridis')

def main():
    print("--- PIPELINE SCRIPT 01a: GSE131907 PREPROCESSING ---")
    
    raw_path = os.path.join(DATA_DIR, 'GSE131907_RAW_merged.h5ad')
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f"Missing raw data: {raw_path}")

    print("1. Loading raw data & running QC...")
    adata = sc.read_h5ad(raw_path)
    adata.X = sp.csr_matrix(adata.X)
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Generate Supp Fig S1: QC Violins
    print("Generating Supplementary Figure S1 (QC Metrics)...")
    adata.obs.rename(columns={'n_genes_by_counts': 'Gene Counts', 'total_counts': 'Total Molecules', 'pct_counts_mt': 'Mitochondrial %'}, inplace=True)
    fig_s1_path = os.path.join(FIG_DIR, 'Supplementary Figures', 'S1', 'Supplementary_Figure_S1_QC.pdf')
    
    # Use return_fig=True to prevent Scanpy from interfering with matplotlib's canvas closing
    fig_s1 = sc.pl.violin(adata, ['Gene Counts', 'Total Molecules', 'Mitochondrial %'], jitter=0.4, multi_panel=True, show=False, return_fig=True)
    fig_s1.savefig(fig_s1_path, bbox_inches='tight')
    plt.close(fig_s1)

    # Revert columns for standard scanpy downstream compatibility
    adata.obs.rename(columns={'Gene Counts': 'n_genes_by_counts', 'Total Molecules': 'total_counts', 'Mitochondrial %': 'pct_counts_mt'}, inplace=True)

    print("2. Applying biological filters...")
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[(adata.obs.n_genes_by_counts < 5000) & (adata.obs.pct_counts_mt < 10), :].copy()
    
    print("3. Normalizing, PCA, and Clustering...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=40, random_state=42)
    sc.tl.umap(adata_hvg, random_state=42)
    # Warning: Reviewers must use the exact same versions of scanpy/leidenalg specified in requirements.txt to reproduce these specific cluster numbers
    sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)
    
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    
    print("4. Annotating Cell Types...")
    cluster_annotations = {
        '1': 'T-Cells', '2': 'Macrophages',
        '3': 'NK / Cytotoxic Cells', '4': 'B-Cells'
    }
    for cluster in ['9', '11', '12', '13', '20', '21', '24', '25', '26']:
        cluster_annotations[cluster] = 'Epithelial cells'
        
    adata.obs['cell_type'] = adata.obs['leiden'].map(lambda x: cluster_annotations.get(x, f"Cluster {x}")).astype('category')
    
    # Generate Supp Fig S2: Marker Dotplot
    print("Generating Supplementary Figure S2 (Marker Dotplot)...")
    marker_genes = ['EPCAM', 'PTPRC', 'CD3D', 'TYROBP', 'NKG7', 'CD79A']
    available_markers = [g for g in marker_genes if g in adata.var_names]
    if available_markers:
        fig_s2_path = os.path.join(FIG_DIR, 'Supplementary Figures', 'S2', 'Supplementary_Figure_S2_Markers.pdf')
        fig_s2 = sc.pl.dotplot(adata, available_markers, groupby='cell_type', standard_scale='var', show=False, return_fig=True)
        fig_s2.savefig(fig_s2_path, bbox_inches='tight')
        plt.close()

    print("5. Mapping Clinical Status & Saving Figure 1...")
    origin_mapping = {
        'nLung': 'Normal', 'nLN': 'Normal',
        'tLung': 'Primary_Sensitive', 'tL/B': 'Primary_Sensitive',
        'mLN': 'Metastatic_Resistant', 'PE': 'Metastatic_Resistant',
        'mBrain': 'Metastatic_Resistant'
    }
    adata.obs['Resistance_Status'] = adata.obs['Sample_Origin'].map(origin_mapping)
    
    # Figure 1: UMAP Exports
    fig1a_path = os.path.join(FIG_DIR, 'Fig 1', 'Figure_1A_CellTypes.pdf')
    fig1a = sc.pl.umap(adata, color='cell_type', palette=sns.color_palette("colorblind").as_hex(), show=False, return_fig=True)
    fig1a.savefig(fig1a_path, bbox_inches='tight')
    plt.close(fig1a)
    
    fig1b_path = os.path.join(FIG_DIR, 'Fig 1', 'Figure_1B_ClinicalOrigin.pdf')
    fig1b = sc.pl.umap(adata, color='Resistance_Status', palette=sns.color_palette("colorblind").as_hex(), show=False, return_fig=True)
    fig1b.savefig(fig1b_path, bbox_inches='tight')
    plt.close(fig1b)

    # Final Export
    final_h5ad = os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad')
    print(f"Saving fully annotated GSE131907 object to: {final_h5ad}")
    adata.write(final_h5ad)
    
    print("--- SCRIPT 01a COMPLETE ---")

if __name__ == "__main__":
    main()
