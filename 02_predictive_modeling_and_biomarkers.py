import scanpy as sc
import pandas as pd
import gseapy as gp
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

def main():
    print("--- Phase 2: Core scRNA-seq Analysis & Modeling ---")
    data_dir = '../data/'
    adata = sc.read_h5ad(data_dir + 'GSE131907_PHASE1_ANNOTATED.h5ad')
    
    print("1. Differential Gene Expression & Artifact Cleaning...")
    tumor_adata = adata[adata.obs['cell_type'] == 'Tumor Cells'].copy()
    
    def group_status(x):
        if x in ['tLung', 'tL/B']: return 'Primary'
        elif x in ['mLN', 'PE', 'mBrain']: return 'Metastatic'
        return 'Other'
    
    tumor_adata.obs['Resistance_Status'] = tumor_adata.obs['Sample_Origin'].map(group_status).astype('category')
    tumor_adata = tumor_adata[tumor_adata.obs['Resistance_Status'].isin(['Primary', 'Metastatic'])].copy()
    
    sc.tl.rank_genes_groups(tumor_adata, groupby='Resistance_Status', reference='Primary', method='wilcoxon')
    df_dge = sc.get.rank_genes_groups_df(tumor_adata, group='Metastatic')
    
    # Clean artifacts for publication
    df_pub = df_dge.rename(columns={'names': 'Gene_Symbol', 'scores': 'Statistical_Score', 'logfoldchanges': 'Log2_Fold_Change', 'pvals_adj': 'Adjusted_P_Value'})
    df_clean = df_pub[(df_pub['Adjusted_P_Value'] < 0.05) & (df_pub['Log2_Fold_Change'] > 0.5) & (df_pub['Log2_Fold_Change'] < 10)]
    df_clean.to_csv(data_dir + 'Supplementary_Data_1_DGE_Cleaned.csv', index=False)
    top_genes = df_clean['Gene_Symbol'].head(150).tolist()

    print("2. Pathway & Master Regulator Analysis...")
    enr_path = gp.enrichr(gene_list=top_genes, gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'], organism='human')
    enr_path.results.to_csv(data_dir + 'Supplementary_Data_2_Pathways.csv', index=False)
    
    enr_tf = gp.enrichr(gene_list=top_genes, gene_sets=['ChEA_2016'], organism='human')
    sig_tfs = enr_tf.results[enr_tf.results['Adjusted P-value'] < 0.05]
    print(f"Discovered {len(sig_tfs)} significant Master Regulators.")
    
    print("3. Random Forest Predictive Modeling...")
    valid_genes = [g for g in top_genes if g in tumor_adata.var_names][:50]
    if 'PCSK1N' not in valid_genes and 'PCSK1N' in tumor_adata.var_names:
        valid_genes[0] = 'PCSK1N'
        
    features = tumor_adata[:, valid_genes].X.toarray() if hasattr(tumor_adata.X, 'toarray') else tumor_adata[:, valid_genes].X
    labels = tumor_adata.obs['Resistance_Status']
    
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)
    rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)
    predictions = rf_model.predict(X_test)
    print("RF Classifier Accuracy Built!")
    
    print("4. Calculating Evolutionary Timeline (Pseudotime)...")
    sc.pp.neighbors(tumor_adata, n_neighbors=15, use_rep='X_pca')
    sc.tl.diffmap(tumor_adata)
    root_id = np.flatnonzero(tumor_adata.obs['Resistance_Status'] == 'Primary')[0]
    tumor_adata.uns['iroot'] = root_id
    sc.tl.dpt(tumor_adata)
    
    tumor_adata.write(data_dir + 'GSE131907_PHASE2_TUMORS.h5ad')
    print("--- Phase 2 Complete ---")

if __name__ == "__main__":
    main()