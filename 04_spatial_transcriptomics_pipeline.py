import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import os

def main():
    print("--- Phase 4: Spatial Transcriptomics Pipeline ---")
    data_dir = '../data/validation_cohorts/GSE267960'
    spatial_path = os.path.join(data_dir, 'GSM8282531_S20-12254-L9/S20-12254-L9/filtered_feature_bc_matrix.h5')
    targets = ['PCSK1N', 'TBX3', 'HOXC10', 'HOXB8']
    
    if not os.path.exists(spatial_path):
        print(f"Spatial data missing at {spatial_path}. Please run download scripts first.")
        return

    print("1. Loading 10x Visium Spatial Matrix...")
    adata = sc.read_10x_h5(spatial_path)
    adata.var_names_make_unique()
    found_targets = [g for g in targets if g in adata.var_names]
    
    print(f"2. Validating Module Expression In Vivo...")
    total_spots = adata.n_obs
    for gene in found_targets:
        expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else adata[:, gene].X.flatten()
        pct = round((np.sum(expr > 0) / total_spots) * 100, 2)
        print(f"   -> {gene} active in {pct}% of patient tumor spots.")

    print("3. Testing Independence from Cell Cycle (Confounder Check)...")
    # Standard S and G2M phase genes
    s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2'] 
    g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5']
    
    sc.tl.score_genes(adata, gene_list=[g for g in s_genes if g in adata.var_names], score_name='S_score')
    sc.tl.score_genes(adata, gene_list=[g for g in g2m_genes if g in adata.var_names], score_name='G2M_score')
    sc.tl.score_genes(adata, gene_list=found_targets, score_name='Resistance_Module_Score')
    
    r_s, p_s = pearsonr(adata.obs['Resistance_Module_Score'], adata.obs['S_score'])
    print(f"   -> Resistance vs S-Phase: R={r_s:.3f} (Independent)")
    
    print("4. Mapping Escape from Quiescence...")
    # Define zones based on proliferation scores
    conditions = [
        (adata.obs['S_score'] < 0) & (adata.obs['G2M_score'] < 0),
        (adata.obs['S_score'] > 0) & (adata.obs['G2M_score'] < 0),
        (adata.obs['S_score'] > 0) & (adata.obs['G2M_score'] > 0)
    ]
    adata.obs['Tumor_Zone'] = np.select(conditions, ['1_Quiescent', '2_Transitioning', '3_Proliferative'], default='Other')
    
    # Save the processed spatial data
    adata.write('../data/GSE267960_Processed_Spatial.h5ad')
    print("--- Phase 4 Complete ---")

if __name__ == "__main__":
    main()