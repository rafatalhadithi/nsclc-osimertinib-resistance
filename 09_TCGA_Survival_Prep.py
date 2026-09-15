import pandas as pd
import numpy as np
import os

# PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))

def main():
    print("--- PIPELINE SCRIPT 09: TCGA SURVIVAL DATA PREPARATION ---")
    
    # 1. Load ML Features
    features_df = pd.read_csv(os.path.join(DATA_DIR, 'Reviewer_Feature_Importances.csv'))
    top_genes = features_df['Gene'].head(10).tolist()
    
    # 2. TCGA Paths
    tcga_dir = os.path.join(DATA_DIR, 'tcga_luad', 'luad_tcga_pan_can_atlas_2018')
    expr_path = os.path.join(tcga_dir, 'data_mrna_seq_v2_rsem.txt')
    clin_path = os.path.join(tcga_dir, 'data_clinical_patient.txt')
    egfr_path = os.path.join(DATA_DIR, 'egfr_status.txt')

    if not os.path.exists(expr_path):
        print(f"[Warning] TCGA dataset not found locally. Skipping extraction.")
        return

    # 3. Clinical & EGFR Filtering
    df_clin = pd.read_csv(clin_path, sep='\t', skiprows=4)
    df_clin = df_clin[['PATIENT_ID', 'PFS_MONTHS', 'PFS_STATUS']].dropna()
    df_clin['Event'] = df_clin['PFS_STATUS'].apply(lambda x: 1 if '1' in str(x) else 0)
    df_clin['Time'] = pd.to_numeric(df_clin['PFS_MONTHS'], errors='coerce')
    
    df_egfr = pd.read_csv(egfr_path, sep='\t')
    egfr_patients = df_egfr[df_egfr.iloc[:, 1] == 1].iloc[:, 0].astype(str).str.split(':').str[-1].str[:12].unique()
    df_clin = df_clin[df_clin['PATIENT_ID'].isin(egfr_patients)]

    # 4. Expression Integration & Z-Scores
    df_expr = pd.read_csv(expr_path, sep='\t')
    matched_genes = [g for g in top_genes if g in df_expr['Hugo_Symbol'].values]
    
    df_expr.set_index('Hugo_Symbol', inplace=True)
    df_t = df_expr.loc[matched_genes].T.reset_index()
    df_t['PATIENT_ID'] = df_t['index'].str[:12]
    
    df = pd.merge(df_clin, df_t, on='PATIENT_ID', how='inner')
    z_scores = df[matched_genes].apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    df['Composite_Z'] = z_scores.mean(axis=1)

    # 5. Quartile Splitting
    q75, q25 = df['Composite_Z'].quantile(0.75), df['Composite_Z'].quantile(0.25)
    high_group = df[df['Composite_Z'] >= q75].assign(Group='High')
    low_group = df[df['Composite_Z'] <= q25].assign(Group='Low')
    
    export_df = pd.concat([high_group, low_group])[['Time', 'Event', 'Group']]
    export_path = os.path.join(DATA_DIR, 'Figure_8_Survival_Data.csv')
    export_df.to_csv(export_path, index=False)
    
    print(f"Data exported for R script: {export_path}")
    print("--- SCRIPT 09 COMPLETE ---")

if __name__ == "__main__":
    main()
