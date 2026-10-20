import os
import sys
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
os.makedirs(DATA_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 02: MAYNARD VALIDATION COHORT PREPROCESSING ---")
    
    meta_path = os.path.join(DATA_DIR, "S01_metacells.csv")
    expr_path = os.path.join(DATA_DIR, "S01_datafinal.csv")
    
    if not os.path.exists(meta_path) or not os.path.exists(expr_path):
        print(f"[ERROR] Missing Maynard raw CSV files in {DATA_DIR}.")
        sys.exit(1)

    df_meta = pd.read_csv(meta_path)
    stage_map = {"naive": "TN", "tn": "TN", "pr": "RD", "rd": "RD", "residual": "RD", "pd": "PD", "progression": "PD"}
    
    stage_col = next((c for c in ["analysis", "treatment_status", "stage"] if c in df_meta.columns), None)
    if stage_col:
        df_meta["Clinical_Stage"] = df_meta[stage_col].astype(str).str.strip().str.lower().map(lambda x: stage_map.get(x, pd.NA))
        
    df_meta_dedup = df_meta.dropna(subset=["cell_id"]).drop_duplicates(subset=["cell_id"]).set_index("cell_id")

    expr_df = pd.read_csv(expr_path, index_col=0)
    if expr_df.shape[0] > expr_df.shape[1] or "SCGB3A2" in expr_df.index: expr_df = expr_df.T 

    adata = sc.AnnData(X=expr_df.values, var=pd.DataFrame(index=expr_df.columns))
    adata.obs_names = [str(x).strip() for x in expr_df.index]
    adata.obs = adata.obs.join(df_meta_dedup, how="left")
    
    adata = adata[adata.obs['Clinical_Stage'].notna()].copy()
    adata.obs['Clinical_Stage'] = adata.obs['Clinical_Stage'].astype('category')
    
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)
    
    adata.write(os.path.join(DATA_DIR, 'maynard_annotated.h5ad'))
    
    # Generate Supp Data 11 (Maynard Manifest)
    patient_col = 'patient_id' if 'patient_id' in adata.obs.columns else adata.obs.columns[0]
    manifest = adata.obs.groupby([patient_col, 'Clinical_Stage']).size().reset_index(name='Retained_Cell_Count')
    manifest.rename(columns={patient_col: 'Patient_ID', 'Clinical_Stage': 'Clinical_Response_State'}, inplace=True)
    manifest['Genotype'] = 'EGFR-mutant'
    manifest['Drug'] = 'Osimertinib'
    manifest['Treatment_Line'] = 'TKI'
    manifest = manifest[['Patient_ID', 'Clinical_Response_State', 'Genotype', 'Drug', 'Treatment_Line', 'Retained_Cell_Count']]
    manifest.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_11.csv'), index=False)

    print("--- SCRIPT 02 COMPLETE ---")

if __name__ == "__main__":
    main()
