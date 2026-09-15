import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# 1. PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
os.makedirs(DATA_DIR, exist_ok=True)

def clean_str(val):
    """Clean string for exact comparison (alphanumeric uppercase only)."""
    if pd.isna(val): return ""
    return "".join(c for c in str(val).strip().upper() if c.isalnum())

def main():
    print("--- PIPELINE SCRIPT 02: MAYNARD VALIDATION COHORT PREPROCESSING ---")
    
    # Files needed in the data directory
    meta_path = os.path.join(DATA_DIR, "S01_metacells.csv")
    expr_path = os.path.join(DATA_DIR, "S01_datafinal.csv")
    
    if not os.path.exists(meta_path) or not os.path.exists(expr_path):
        print(f"[ERROR] Missing required Maynard CSV files in {DATA_DIR}.")
        print("Please ensure 'S01_metacells.csv' and 'S01_datafinal.csv' are present.")
        sys.exit(1)

    print("1. Loading metadata...")
    df_meta = pd.read_csv(meta_path)
    
    # Standardize clinical stage
    stage_map = {
        "naive": "TN", "naïve": "TN", "na•ve": "TN", "tn": "TN",
        "grouped_pr": "RD", "pr": "RD", "rd": "RD", "residual": "RD",
        "grouped_pd": "PD", "pd": "PD", "progression": "PD", "post-osimertinib": "RD"
    }
    
    # Locate analysis/stage column
    stage_col = next((c for c in ["analysis", "treatment_status", "stage"] if c in df_meta.columns), None)
    if stage_col:
        raw_stage = df_meta[stage_col].astype(str).str.strip().str.lower()
        df_meta["Clinical_Stage"] = raw_stage.map(lambda x: stage_map.get(x, np.nan))
    else:
        df_meta["Clinical_Stage"] = np.nan

    # Deduplicate metadata to ensure 1-to-1 match with expression matrix
    df_meta_dedup = df_meta.dropna(subset=["cell_id"]).drop_duplicates(subset=["cell_id"])
    df_meta_dedup.set_index("cell_id", inplace=True)

    print("2. Loading and orienting expression matrix...")
    expr_df = pd.read_csv(expr_path, index_col=0)
    
    # Transpose if genes are detected as rows
    if expr_df.shape[0] > expr_df.shape[1] or "SCGB3A2" in expr_df.index:
        expr_df = expr_df.T 

    print("3. Constructing AnnData Object...")
    adata = sc.AnnData(X=expr_df.values, var=pd.DataFrame(index=expr_df.columns))
    adata.obs_names = [str(x).strip() for x in expr_df.index]
    
    # Merge metadata
    adata.obs = adata.obs.join(df_meta_dedup, how="left")
    
    # Filter to valid clinical stages
    adata = adata[adata.obs['Clinical_Stage'].notna()].copy()
    adata.obs['Clinical_Stage'] = adata.obs['Clinical_Stage'].astype('category')
    
    # Ensure tumor isolation (If general_annotation column exists)
    if 'general_annotation' in adata.obs.columns:
        tumor_mask = adata.obs['general_annotation'].astype(str).str.contains('Epithelial|Cancer|Tumor', case=False, na=False)
        if tumor_mask.sum() > 0:
            adata = adata[tumor_mask].copy()

    print("4. Basic Quality Control and Normalization...")
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)
    # Note: Log1p is applied downstream dynamically where needed, storing raw counts here
    
    print(f"Final Validation Dataset: {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs['Clinical_Stage'].value_counts())

    # Final Export
    final_h5ad = os.path.join(DATA_DIR, 'maynard_annotated.h5ad')
    print(f"5. Saving Maynard annotated object to: {final_h5ad}")
    adata.write(final_h5ad)
    print("--- SCRIPT 02 COMPLETE ---")

if __name__ == "__main__":
    main()
