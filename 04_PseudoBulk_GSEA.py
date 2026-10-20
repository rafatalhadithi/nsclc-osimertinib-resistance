import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 2'))
os.makedirs(FIG_DIR, exist_ok=True)

def run_pseudobulk_gsea(adata, patient_col, group_col, ref_group, target_group, output_prefix):
    pdata = dc.get_pseudobulk(adata, sample_col=patient_col, groups_col=group_col, mode='sum', min_cells=5, min_counts=500, skip_checks=True)
    pdata.X = np.clip(pdata.X, a_min=0, a_max=None)
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)
    
    pdata.obs[group_col] = pdata.obs[group_col].astype('category')
    sc.tl.rank_genes_groups(pdata, groupby=group_col, reference=ref_group, method='wilcoxon')
    deg_df = sc.get.rank_genes_groups_df(pdata, group=target_group).dropna(subset=['logfoldchanges', 'names'])
    
    if output_prefix == "Discovery":
        sig_df = deg_df[(deg_df['logfoldchanges'] > 0.25) & (deg_df['pvals_adj'] < 0.1)].sort_values('logfoldchanges', ascending=False)
        sig_df.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_1.csv'), index=False)
    
    deg_df['rank_metric'] = deg_df['logfoldchanges']
    deg_df = deg_df.sort_values('rank_metric', ascending=False)
    rnk_path = os.path.join(DATA_DIR, f'{output_prefix}_Ranked.rnk')
    deg_df[['names', 'rank_metric']].to_csv(rnk_path, sep='\t', index=False, header=False)
    
    gsea_res = gp.prerank(rnk=rnk_path, gene_sets='MSigDB_Hallmark_2020', threads=4, min_size=5, max_size=500, permutation_num=1000, seed=42)
    out_gsea = os.path.join(DATA_DIR, 'Supplementary_Data_2.csv' if output_prefix == "Discovery" else 'Supplementary_Data_10.csv')
    gsea_res.res2d.to_csv(out_gsea)
    return gsea_res.res2d

def main():
    print("--- PIPELINE SCRIPT 04: PSEUDO-BULK & GSEA (FIGURE 2) ---")
    adata_disc = sc.read_h5ad(os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad'))
    if adata_disc.raw is not None: adata_disc = adata_disc.raw.to_adata()
    
    cell_type_col = next((c for c in adata_disc.obs.columns if 'cell' in c.lower()), 'cell_type')
    origin_col = next((c for c in adata_disc.obs.columns if 'origin' in c.lower() or 'site' in c.lower()), 'Sample_Origin')
    sample_col = next((c for c in adata_disc.obs.columns if 'patient' in c.lower() or 'sample' in c.lower()), 'Sample')
    
    adata_disc_tumor = adata_disc[adata_disc.obs[cell_type_col].astype(str).str.contains('Epithelial|Tumor', case=False, na=False)].copy()
    adata_disc_tumor.obs['Status'] = adata_disc_tumor.obs[origin_col].map({'tLung': 'Primary', 'tL/B': 'Primary', 'PE': 'Metastatic', 'mLN': 'Metastatic'})
    adata_disc_tumor = adata_disc_tumor[adata_disc_tumor.obs['Status'].notna()].copy()
    
    disc_summary = run_pseudobulk_gsea(adata_disc_tumor, sample_col, 'Status', 'Primary', 'Metastatic', 'Discovery')
    
    adata_val = sc.read_h5ad(os.path.join(DATA_DIR, 'maynard_annotated.h5ad'))
    if adata_val.raw is not None: adata_val = adata_val.raw.to_adata()
    
    val_sub = adata_val[adata_val.obs['Clinical_Stage'].isin(['TN', 'PD'])].copy()
    patient_col_val = next((col for col in ['patient_id', 'sample_id', 'Patient'] if col in val_sub.obs.columns), val_sub.obs.columns[0])
    
    val_summary = run_pseudobulk_gsea(val_sub, patient_col_val, 'Clinical_Stage', 'TN', 'PD', 'Maynard')
    
    disc_summary['Term'] = disc_summary['Term'].str.replace('Hallmark ', '', regex=False).str.replace('_', ' ', regex=False).replace('heme Metabolism', 'Heme Metabolism')
    val_summary['Term'] = val_summary['Term'].str.replace('Hallmark ', '', regex=False).str.replace('_', ' ', regex=False).replace('heme Metabolism', 'Heme Metabolism')
    
    target_pathways = ['Oxidative Phosphorylation', 'Glycolysis', 'Myc Targets V1', 'Interferon Alpha Response', 'E2F Targets', 'Epithelial Mesenchymal Transition', 'Angiogenesis', 'Pancreas Beta Cells', 'KRAS Signaling Dn', 'Heme Metabolism']
    df_disc_sub = disc_summary[disc_summary['Term'].isin(target_pathways)][['Term', 'NES']].assign(Cohort='Discovery (GSE131907)')
    df_val_sub = val_summary[val_summary['Term'].isin(target_pathways)][['Term', 'NES']].assign(Cohort='Validation (Maynard PD vs TN)')
    
    plot_df = pd.concat([df_disc_sub, df_val_sub], ignore_index=True)
    plot_df['Term'] = pd.Categorical(plot_df['Term'], categories=target_pathways, ordered=True)
    
    plt.figure(figsize=(11, 7))
    sns.barplot(data=plot_df.sort_values('Term'), x='NES', y='Term', hue='Cohort', palette={'Discovery (GSE131907)': '#4682B4', 'Validation (Maynard PD vs TN)': '#DC143C'}, edgecolor='black')
    plt.axvline(x=0, color='black', linewidth=1.5); plt.ylabel('')
    plt.title('Figure 2: Cross-Cohort Validation of Hallmark Resistance Pathways', fontweight='bold', pad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_2_Pathways_Cross_Validated.pdf'), format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    print("--- SCRIPT 04 COMPLETE ---")

if __name__ == "__main__":
    main()
