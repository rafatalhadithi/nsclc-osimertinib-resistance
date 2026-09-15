import os
import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 1. PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures/Fig 2'))
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

def run_pseudobulk_gsea(adata, patient_col, group_col, ref_group, target_group, output_prefix):
    """Aggregates to pseudobulk, calculates DEGs, and runs GSEA."""
    print(f"Aggregating {output_prefix} by {patient_col} and {group_col}...")
    
    # Decoupler Pseudo-bulking
    pdata = dc.pp.pseudobulk(adata, sample_col=patient_col, groups_col=group_col, mode='sum', skip_checks=True)
    dc.pp.filter_samples(pdata, min_cells=5, min_counts=500)
    
    # Normalization & DE
    pdata.X = np.clip(pdata.X, a_min=0, a_max=None)
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)
    
    sc.tl.rank_genes_groups(pdata, groupby=group_col, reference=ref_group, method='wilcoxon')
    deg_df = sc.get.rank_genes_groups_df(pdata, group=target_group).dropna(subset=['logfoldchanges', 'names'])
    
    # Save DGE for manuscript (Supplementary Data 1)
    if output_prefix == "Discovery":
        sig_df = deg_df[(deg_df['logfoldchanges'] > 0.25) & (deg_df['pvals_adj'] < 0.1)].sort_values('logfoldchanges', ascending=False)
        sig_df.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_1_DGE_Cleaned.csv'), index=False)
    
    # Continuous Log2FC ranking to prevent GSEA ties
    deg_df['rank_metric'] = deg_df['logfoldchanges']
    deg_df = deg_df.sort_values('rank_metric', ascending=False)
    
    rnk_path = os.path.join(DATA_DIR, f'{output_prefix}_Ranked.rnk')
    deg_df[['names', 'rank_metric']].to_csv(rnk_path, sep='\t', index=False, header=False)
    
    print(f"Running GSEA for {output_prefix}...")
    gsea_res = gp.prerank(rnk=rnk_path, gene_sets='MSigDB_Hallmark_2020', threads=4, min_size=5, max_size=500, permutation_num=1000, seed=42)
    
    summary_df = gsea_res.res2d
    summary_df.to_csv(os.path.join(DATA_DIR, f'{output_prefix}_GSEA_Summary.csv'))
    return summary_df

def main():
    print("--- PIPELINE SCRIPT 04: PSEUDO-BULK & GSEA (FIGURE 2) ---")
    
    # ---------------------------------------------------------
    # PART A: DISCOVERY COHORT (GSE131907)
    # ---------------------------------------------------------
    print("\n1. Processing Discovery Cohort...")
    adata_disc = sc.read_h5ad(os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad'))
    if adata_disc.raw is not None: adata_disc = adata_disc.raw.to_adata()
    
    tumor_mask = adata_disc.obs['cell_type'].astype(str).str.contains('Epithelial', case=False, na=False)
    adata_disc_tumor = adata_disc[tumor_mask].copy()
    
    disc_summary = run_pseudobulk_gsea(
        adata_disc_tumor, patient_col='Sample', group_col='Resistance_Status', 
        ref_group='Primary_Sensitive', target_group='Metastatic_Resistant', output_prefix='Discovery'
    )
    
    # ---------------------------------------------------------
    # PART B: VALIDATION COHORT (MAYNARD)
    # ---------------------------------------------------------
    print("\n2. Processing Validation Cohort...")
    adata_val = sc.read_h5ad(os.path.join(DATA_DIR, 'maynard_annotated.h5ad'))
    if adata_val.raw is not None: adata_val = adata_val.raw.to_adata()
    
    val_summary = run_pseudobulk_gsea(
        adata_val, patient_col='sample_id', group_col='Clinical_Stage', 
        ref_group='TN', target_group='PD', output_prefix='Maynard_PD_vs_TN'
    )
    
    # ---------------------------------------------------------
    # PART C: PLOTTING FIGURE 2 (CROSS-COHORT GSEA)
    # ---------------------------------------------------------
    print("\n3. Generating Figure 2 (Cross-Cohort GSEA)...")
    target_pathways = [
        'Heme Metabolism', 'Pancreas Beta Cells', 'Cholesterol Homeostasis', 
        'Peroxisome', 'KRAS Signaling Dn', 'Interferon Gamma Response', 
        'Inflammatory Response', 'UV Response Dn', 'p53 Pathway', 'TNF-alpha Signaling via NF-kB'
    ]
    
    # Clean terms
    disc_summary['Term'] = disc_summary['Term'].str.replace('Hallmark ', '', regex=False).str.replace('_', ' ', regex=False).replace('Pperoxisome', 'Peroxisome').replace('heme Metabolism', 'Heme Metabolism')
    val_summary['Term'] = val_summary['Term'].str.replace('Hallmark ', '', regex=False).str.replace('_', ' ', regex=False).replace('Pperoxisome', 'Peroxisome').replace('heme Metabolism', 'Heme Metabolism')
    
    df_disc_sub = disc_summary[disc_summary['Term'].isin(target_pathways)][['Term', 'NES']].assign(Cohort='Discovery (GSE131907)')
    df_val_sub = val_summary[val_summary['Term'].isin(target_pathways)][['Term', 'NES']].assign(Cohort='Validation (Maynard PD vs TN)')
    
    plot_df = pd.concat([df_disc_sub, df_val_sub], ignore_index=True)
    plot_df['Term'] = pd.Categorical(plot_df['Term'], categories=target_pathways, ordered=True)
    plot_df = plot_df.sort_values('Term')
    
    plt.figure(figsize=(11, 7))
    sns.set_theme(style="ticks", context="paper")
    ax = sns.barplot(data=plot_df, x='NES', y='Term', hue='Cohort', palette={'Discovery (GSE131907)': '#4682B4', 'Validation (Maynard PD vs TN)': '#DC143C'}, edgecolor='black', linewidth=0.8)
    
    plt.title('Figure 2: Cross-Cohort Validation of Hallmark Resistance Pathways', fontsize=15, fontweight='bold', pad=15)
    plt.xlabel('Normalized Enrichment Score (NES)', fontsize=13, fontweight='bold')
    plt.ylabel('')
    plt.axvline(x=0, color='black', linewidth=1.5)
    plt.legend(title='Clinical Cohort', loc='lower right', framealpha=0.9)
    plt.tight_layout()
    
    fig2_pdf = os.path.join(FIG_DIR, 'Figure_2_Pathways_Cross_Validated.pdf')
    plt.savefig(fig2_pdf, format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"   -> Saved: {fig2_pdf}")
    print("--- SCRIPT 04 COMPLETE ---")

if __name__ == "__main__":
    main()
