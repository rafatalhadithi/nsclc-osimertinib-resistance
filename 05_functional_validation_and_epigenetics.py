import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

def main():
    print("--- Phase 5: Functional Validation & Epigenetics ---")
    data_dir = '../data/'
    
    print("1. Mining Real CRISPR-Cas9 Synthetic Lethality Data...")
    crispr_file = f"{data_dir}Pfeifer_et_al_Supplemental_Data/MAGECK_files/PC9_Osimertinib_CRISPRn.txt"
    cascade_genes = ['PCSK1N', 'TBX3', 'HOXC10', 'HOXB8']
    
    try:
        df_crispr = pd.read_csv(crispr_file, sep='\t')
        df_crispr.columns = [col.lower() for col in df_crispr.columns]
        gene_col = [c for c in df_crispr.columns if c in ['id', 'gene']][0]
        lfc_col = [c for c in df_crispr.columns if 'lfc' in c][0]
        
        df_crispr[gene_col] = df_crispr[gene_col].astype(str).str.upper()
        res = df_crispr[df_crispr[gene_col].isin(cascade_genes)][[gene_col, lfc_col]]
        res.to_csv(f"{data_dir}CRISPR_Functional_Validation.csv", index=False)
        print("   -> CRISPR Metrics Saved.")
    except Exception as e:
        print(f"   -> Could not load CRISPR data: {e}")

    print("2. Epigenetic Locus Discovery (EZH2-HOXC Seesaw)...")
    spatial_path = f"{data_dir}validation_cohorts/GSE267960/GSM8282531_S20-12254-L9/S20-12254-L9/filtered_feature_bc_matrix.h5"
    if os.path.exists(spatial_path):
        adata = sc.read_10x_h5(spatial_path)
        adata.var_names_make_unique()
        
        # Test spatial inverse correlation
        if 'EZH2' in adata.var_names and 'HOXC9' in adata.var_names:
            ezh2 = adata[:, 'EZH2'].X.toarray().flatten() if hasattr(adata[:, 'EZH2'].X, "toarray") else adata[:, 'EZH2'].X.flatten()
            hoxc9 = adata[:, 'HOXC9'].X.toarray().flatten() if hasattr(adata[:, 'HOXC9'].X, "toarray") else adata[:, 'HOXC9'].X.flatten()
            
            spot_corr, spot_pval = spearmanr(ezh2, hoxc9)
            print(f"   -> True Spatial Correlation (EZH2 vs HOXC9): R = {spot_corr:.3f}")
            
            print("3. Simulating Bulk Masking (Simpson's Paradox)...")
            n_pseudo_patients = 150
            pseudo_ezh2, pseudo_hox = [], []
            for _ in range(n_pseudo_patients):
                idx = np.random.choice(len(ezh2), size=50)
                pseudo_ezh2.append(np.mean(ezh2[idx]))
                pseudo_hox.append(np.mean(hoxc9[idx]))
                
            bulk_corr, bulk_pval = spearmanr(pseudo_ezh2, pseudo_hox)
            print(f"   -> Pseudo-Bulk Correlation (Masked): R = {bulk_corr:.3f}")
            
            pd.DataFrame({
                'Resolution': ['Spatial Micro-Environment', 'Simulated Bulk'],
                'Spearman_R': [spot_corr, bulk_corr]
            }).to_csv(f"{data_dir}Pseudo_Bulk_Masking_Simulation.csv", index=False)

    print("--- Phase 5 Complete ---")

if __name__ == "__main__":
    main()