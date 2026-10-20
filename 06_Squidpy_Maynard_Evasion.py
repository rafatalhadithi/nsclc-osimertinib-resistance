import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 4'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 06: SQUIDPY IMMUNE EVASION (MAYNARD SMART-SEQ2) ---")
    adata = sc.read_h5ad(os.path.join(DATA_DIR, 'maynard_annotated.h5ad'))
    
    stage_col = next((c for c in ['Clinical_Stage', 'stage_clean', 'analysis', 'treatment_status'] if c in adata.obs.columns), None)
    stage_map = {'tn': 'TN', 'naive': 'TN', 'pd': 'PD', 'progression': 'PD', 'progressive disease': 'PD'}
    adata.obs['Mapped_Stage'] = adata.obs[stage_col].astype(str).str.lower().str.strip().map(lambda x: stage_map.get(x, 'Other'))
    adata = adata[adata.obs['Mapped_Stage'].isin(['TN', 'PD'])].copy()
    
    cell_col = next((c for c in ['general_annotation', 'cell_type', 'Cell_type'] if c in adata.obs.columns), None)
    
    def group_status(row):
        ctype = str(row[cell_col]).lower().replace('_', ' ')
        stage = str(row['Mapped_Stage'])
        if any(k in ctype for k in ['epi', 'tumor', 'cancer']): return f'Tumor ({stage})'
        if any(k in ctype for k in ['macrophage', 'myeloid', 'monocyte', 'dc']): return 'Myeloid/Macrophages'
        if any(k in ctype for k in ['t cell', 'nk', 'lymphocyte']): return 'T/NK Cells'
        if 'immune' in ctype: return 'Immune Cells'
        return 'Other'
        
    adata.obs['Interaction_Group'] = adata.obs.apply(group_status, axis=1)
    adata_clean = adata[adata.obs['Interaction_Group'] != 'Other'].copy()
    adata_clean.obs['Interaction_Group'] = adata_clean.obs['Interaction_Group'].astype('category')
    
    sc.pp.normalize_total(adata_clean, target_sum=1e4)
    sc.pp.log1p(adata_clean)

    res = sq.gr.ligrec(adata_clean, cluster_key='Interaction_Group', n_perms=1000, threshold=0.0, copy=True, seed=42, use_raw=False)
    means_df, pvals_df = res['means'], res['pvalues']
    
    target_interactions = [('THBS1', 'CD47'), ('CD274', 'PDCD1'), ('LGALS9', 'HAVCR2'), ('PVR', 'TIGIT'), ('NECTIN2', 'TIGIT'), ('MIF', 'CD74'), ('CD86', 'CTLA4')]
    selected_cols = [(src, tgt) for src in ['Tumor (TN)', 'Tumor (PD)'] for tgt in ['Myeloid/Macrophages', 'T/NK Cells', 'Immune Cells'] if (src, tgt) in pvals_df.columns]
    
    plot_data = []
    for lig, rec in target_interactions:
        actual_pair = (lig, rec) if (lig, rec) in means_df.index else (rec, lig) if (rec, lig) in means_df.index else None
        if actual_pair:
            for col in selected_cols:
                pval, mean_exp = pvals_df.loc[actual_pair, col], means_df.loc[actual_pair, col]
                if pd.notna(mean_exp) and mean_exp > 0 and pd.notna(pval) and pval <= 0.05:
                    plot_data.append({'Pair': f"{lig} - {rec}", 'Interaction': f"{col[0]} -> \n{col[1]}", 'p_val': float(pval), 'neg_log_p': -np.log10(float(pval) + 1e-4), 'Mean_Expression': float(mean_exp)})
                    
    df_plot = pd.DataFrame(plot_data)
    if not df_plot.empty:
        df_plot = df_plot.sort_values(by=['Pair', 'Interaction'])
        fig, ax = plt.subplots(figsize=(8, 5))
        scatter = ax.scatter(x=df_plot['Interaction'], y=df_plot['Pair'], s=df_plot['neg_log_p'] * 150 + 50, c=df_plot['Mean_Expression'], cmap='plasma', edgecolors='black', linewidth=1.2, zorder=3)
        for _, row in df_plot.iterrows():
            sig = "***" if row['p_val'] <= 0.001 else "**" if row['p_val'] <= 0.01 else "*" if row['p_val'] <= 0.05 else ""
            if sig: ax.annotate(sig, (row['Interaction'], row['Pair']), xytext=(15, -3), textcoords="offset points", ha='left', va='center', fontsize=12, fontweight='bold', zorder=4)

        ax.grid(True, linestyle='--', alpha=0.4, zorder=1); ax.set_axisbelow(True)
        plt.xticks(rotation=45, ha='right', fontsize=10); plt.margins(0.2)
        plt.title('Tumor-Immune Checkpoint Interactions\n(Smart-seq2 Cohort, N=1,000 Permutations)', pad=20, fontweight='bold')
        plt.colorbar(scatter, ax=ax, pad=0.04).set_label('Mean Interaction Strength', rotation=270, labelpad=18)
        
        plt.savefig(os.path.join(FIG_DIR, 'Figure_4_Immune_Evasion.pdf'), dpi=600, transparent=True, bbox_inches='tight')
        plt.close()
    print("--- SCRIPT 06 COMPLETE ---")

if __name__ == "__main__":
    main()
