import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 3'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 05: PAGA TRAJECTORY (COMPACT HIGH-RES NETWORK) ---")
    adata = sc.read_h5ad(os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad'))
    
    cell_col = next((c for c in adata.obs.columns if 'cell' in c.lower()), 'cell_type')
    origin_col = next((c for c in adata.obs.columns if 'origin' in c.lower() or 'site' in c.lower()), 'Sample_Origin')
    
    adata_tumor = adata[adata.obs[cell_col].astype(str).str.contains('Epithelial|Tumor|Cancer', case=False, na=False)].copy()
    
    origin_mapping = {'nLung': 'Normal Lung', 'tLung': 'Primary Tumor (Lung)', 'tL/B': 'Primary Tumor (Bronchus)', 'PE': 'Metastatic Fluid', 'mLN': 'Metastatic Lymph Node', 'mBrain': 'Brain Metastasis'}
    adata_tumor = adata_tumor[adata_tumor.obs[origin_col].isin(origin_mapping.keys())].copy()
    adata_tumor.obs['PAGA_Group'] = adata_tumor.obs[origin_col].map(origin_mapping).astype('category')
    
    color_dict = {'Normal Lung':'#ff9933', 'Primary Tumor (Lung)':'#ffff00', 'Primary Tumor (Bronchus)':'#ffff00', 'Metastatic Fluid':'#800080', 'Metastatic Lymph Node':'#d62728', 'Brain Metastasis':'#000080'}
    
    if 'X_pca' not in adata_tumor.obsm: sc.pp.pca(adata_tumor, svd_solver='arpack', random_state=42)
    sc.pp.neighbors(adata_tumor, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.paga(adata_tumor, groups='PAGA_Group')
    
    adj = adata_tumor.uns['paga']['connectivities'].todense()
    groups = adata_tumor.obs['PAGA_Group'].cat.categories
    cell_counts = adata_tumor.obs['PAGA_Group'].value_counts()[groups].values
    
    G = nx.Graph()
    for i in range(len(groups)): G.add_node(i)
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            if adj[i, j] > 0.05: G.add_edge(i, j, weight=adj[i, j] * 5) 
                
    pos = nx.spring_layout(G, weight='weight', k=0.6, seed=42)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('off') 
    
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            if adj[i, j] > 0.05: ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]], color='#d3d3d3', linewidth=adj[i, j] * 8, zorder=1)
                
    for i, group in enumerate(groups):
        ax.scatter(pos[i][0], pos[i][1], s=cell_counts[i] * 0.2, c=color_dict[group], edgecolors='black', linewidths=1.5, zorder=2)
        ax.text(pos[i][0], pos[i][1] + 0.08, group.replace(' (Lung)', '').replace(' (Bronchus)', ''), ha='center', va='bottom', fontsize=11, fontweight='bold', zorder=3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_3C_PAGA_Network.pdf'), dpi=600, transparent=True, bbox_inches='tight')
    plt.close()
    
    pd.DataFrame(adj, index=groups, columns=groups).to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_3.csv'))
    print("--- SCRIPT 05 COMPLETE ---")

if __name__ == "__main__":
    main()
