import scanpy as sc
import pandas as pd
import scipy.sparse as sp

def main():
    print("--- Phase 1: Preprocessing and Annotation ---")
    data_dir = '../data/' 
    
    print("1. Loading raw data & running QC...")
    adata = sc.read_h5ad(data_dir + 'GSE131907_RAW_merged.h5ad')
    adata.X = sp.csr_matrix(adata.X)
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    print("2. Applying biological filters...")
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    
    print("3. Normalizing, PCA, and Clustering...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    
    # Slice to HVGs for memory-safe scaling
    adata_hvg = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=0.5)
    
    # Sync clustering back to the main object
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    
    print("4. Finding Markers & Annotating Cell Types...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(data_dir + 'marker_genes_complete.csv', index=False)
    
    cluster_annotations = {
        '1': 'T-Cells (CD3D+)', '2': 'Macrophages (TYROBP+)',
        '3': 'NK / Cytotoxic Cells (NKG7+)', '4': 'B-Cells (CD79A+)'
    }
    
    # Include manually curated Tumor Clusters
    tumor_clusters = ['9', '11', '12', '13', '20', '21', '24', '25', '26']
    for cluster in tumor_clusters:
        cluster_annotations[cluster] = 'Tumor Cells'
        
    adata.obs['cell_type'] = adata.obs['leiden'].map(lambda x: cluster_annotations.get(x, f"Cluster {x}"))
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    
    print("5. Mapping Clinical Status...")
    origin_mapping = {
        'nLung': '1. Normal Lung', 'nLN': '2. Normal Lymph Node',
        'tLung': '3. Primary Tumor', 'tL/B': '3. Primary Tumor',
        'mLN': '4. Metastatic Lymph Node', 'PE': '5. Metastatic Fluid',
        'mBrain': '6. Brain Metastasis'
    }
    adata.obs['Clinical_Status'] = adata.obs['Sample_Origin'].map(origin_mapping)
    
    print("Saving final Phase 1 dataset...")
    adata.write(data_dir + 'GSE131907_PHASE1_ANNOTATED.h5ad')
    print("--- Phase 1 Complete ---")

if __name__ == "__main__":
    main()