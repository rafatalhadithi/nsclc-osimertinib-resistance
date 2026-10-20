import os
import sys
import pandas as pd
import scanpy as sc
import matplotlib.subplots as plt
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 6'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 07: MAYNARD VIOLIN PLOTS (FIGURE 6) ---")
    adata = sc.read_h5ad(os.path.join(DATA_DIR, 'maynard_annotated.h5ad'))
    if adata.raw is not None: adata = adata.raw.to_adata()
        
    stage_col = next((c for c in ['Clinical_Stage', 'stage_clean', 'analysis', 'treatment_status'] if c in adata.obs.columns), None)
    stage_map = {'tn': 'TN', 'naive': 'TN', 'rd': 'RD', 'residual': 'RD', 'pr': 'RD', 'pd': 'PD', 'progression': 'PD', 'progressive disease': 'PD'}
    adata.obs['Mapped_Stage'] = adata.obs[stage_col].astype(str).str.lower().str.strip().map(lambda x: stage_map.get(x, 'Other'))
    
    adata = adata[adata.obs['Mapped_Stage'].isin(['TN', 'RD', 'PD'])].copy()
    adata.obs['Mapped_Stage'] = pd.Categorical(adata.obs['Mapped_Stage'], categories=['TN', 'RD', 'PD'], ordered=True)
    
    target_genes = ['SCGB3A2', 'MUC1', 'MUC5B', 'EZH2']
    found_genes = [g for g in target_genes if g in adata.var_names]
    
    df_plot = pd.DataFrame(index=adata.obs.index)
    df_plot['Clinical Stage'] = adata.obs['Mapped_Stage']
    
    for gene in found_genes:
        gene_data = adata[:, gene].X
        df_plot[gene] = gene_data.toarray().flatten() if hasattr(gene_data, "toarray") else gene_data.flatten()
            
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    palette = {'TN': '#4c72b0', 'RD': '#dd8452', 'PD': '#c44e52'}
    
    for i, gene in enumerate(target_genes):
        if gene in found_genes:
            sns.violinplot(data=df_plot, x='Clinical Stage', y=gene, ax=axes[i], palette=palette, scale='width', inner=None, linewidth=1.2)
            axes[i].set_title(f"{gene} Expression", fontweight='bold', fontsize=12)
            axes[i].set_ylabel("Log(Normalized Expression)"); axes[i].set_xlabel("Clinical Stage")
            axes[i].grid(True, axis='y', linestyle='--', alpha=0.7)
            axes[i].spines['top'].set_visible(False); axes[i].spines['right'].set_visible(False)
            
    plt.tight_layout(pad=3.0)
    plt.savefig(os.path.join(FIG_DIR, 'Figure_6_Maynard_Violins.pdf'), dpi=600, bbox_inches='tight')
    plt.close()
    
    df_plot.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_9.csv'))
    print("--- SCRIPT 07 COMPLETE ---")

if __name__ == "__main__":
    main()
