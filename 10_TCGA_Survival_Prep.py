import pandas as pd
import numpy as np
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
TCGA_DIR = os.path.join(DATA_DIR, 'tcga_luad', 'luad_tcga_pan_can_atlas_2018')

def main():
    print("--- PIPELINE SCRIPT 10: TCGA SURVIVAL DATA PREPARATION ---")
    feature_file = os.path.join(DATA_DIR, 'Supplementary_Data_5.csv')
    if not os.path.exists(feature_file): sys.exit(1)

    top_genes = pd.read_csv(feature_file).iloc[:, 0].head(10).tolist()
    
    expr_path, clin_path, egfr_path = os.path.join(TCGA_DIR, 'data_mrna_seq_v2_rsem.txt'), os.path.join(TCGA_DIR, 'data_clinical_patient.txt'), os.path.join(DATA_DIR, 'egfr_status.txt')
    if not os.path.exists(expr_path): sys.exit(1)

    df_clin = pd.read_csv(clin_path, sep='\t', skiprows=4)[['PATIENT_ID', 'OS_MONTHS', 'OS_STATUS', 'PFS_MONTHS', 'PFS_STATUS']].dropna()
    df_clin['OS_Event'] = df_clin['OS_STATUS'].apply(lambda x: 1 if '1' in str(x) else 0)
    df_clin['PFS_Event'] = df_clin['PFS_STATUS'].apply(lambda x: 1 if '1' in str(x) else 0)
    df_clin['OS_Time'] = pd.to_numeric(df_clin['OS_MONTHS'], errors='coerce')
    df_clin['PFS_Time'] = pd.to_numeric(df_clin['PFS_MONTHS'], errors='coerce')
    
    df_egfr = pd.read_csv(egfr_path, sep='\t')
    egfr_col = [col for col in df_egfr.columns if 'EGFR' in col.upper()][0]
    df_egfr[egfr_col] = pd.to_numeric(df_egfr[egfr_col], errors='coerce')
    df_clin = df_clin[df_clin['PATIENT_ID'].isin(df_egfr[df_egfr[egfr_col] == 1].iloc[:, 0].astype(str).str.split(':').str[-1].str[:12].unique())]

    df_expr = pd.read_csv(expr_path, sep='\t').set_index('Hugo_Symbol')
    matched_genes = [g for g in top_genes if g in df_expr.index]
    df_t = df_expr.loc[matched_genes].T.reset_index()
    df_t['PATIENT_ID'] = df_t['index'].str[:12]
    
    df = pd.merge(df_clin, df_t, on='PATIENT_ID', how='inner').dropna(subset=matched_genes)
    for g in matched_genes: df[g] = pd.to_numeric(df[g], errors='coerce')
    
    df['Composite_Z'] = df[matched_genes].apply(lambda x: (x - x.mean()) / x.std(), axis=0).mean(axis=1)

    # PFS 
    q75, q25 = df['Composite_Z'].quantile(0.75), df['Composite_Z'].quantile(0.25)
    export_df_pfs = pd.concat([df[df['Composite_Z'] >= q75].assign(Group='High'), df[df['Composite_Z'] <= q25].assign(Group='Low')])[['PFS_Time', 'PFS_Event', 'Group']].rename(columns={'PFS_Time': 'Time', 'PFS_Event': 'Event'})
    export_df_pfs.to_csv(os.path.join(DATA_DIR, 'Figure_9_Survival_Data.csv'), index=False)
    
    # OS
    df['Group'] = np.where(df['Composite_Z'] >= df['Composite_Z'].median(), 'High', 'Low')
    df[['OS_Time', 'OS_Event', 'Group']].rename(columns={'OS_Time': 'Time', 'OS_Event': 'Event'}).to_csv(os.path.join(DATA_DIR, 'Figure_7B_Survival_Data.csv'), index=False)

    print("--- SCRIPT 10 COMPLETE ---")

if __name__ == "__main__":
    main()
