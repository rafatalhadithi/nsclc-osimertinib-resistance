import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
import warnings
warnings.filterwarnings('ignore')

def main():
    print("--- Phase 6: Advanced Rigor & Confounder Ablation ---")
    data_dir = '../data/'
    file_path = f"{data_dir}GSE131907_FINAL_ANNOTATED.h5ad"
    
    print("1. Loading Single-Cell Data for Ablation Testing...")
    adata = sc.read_h5ad(file_path)

    # Isolate Tumor Cells and remove specific confounding cohorts (e.g., mBrain)
    tumor_cells = adata[adata.obs['Cell_type'] == 'Epithelial cells'].copy()
    valid_origins = ['tLung', 'tL/B', 'PE', 'mLN']
    tumor_cells = tumor_cells[tumor_cells.obs['Sample_Origin'].isin(valid_origins)].copy()

    def assign_status(origin):
        if origin in ['tLung', 'tL/B']: return 'Primary_Sensitive'
        elif origin in ['PE', 'mLN']: return 'Metastatic_Resistant'
        return 'Unknown'
    
    tumor_cells.obs['Resistance_Status'] = tumor_cells.obs['Sample_Origin'].apply(assign_status)

    print("2. Running Stratified 5-Fold Cross Validation...")
    if 'highly_variable' not in tumor_cells.var.columns:
        sc.pp.highly_variable_genes(tumor_cells, n_top_genes=50)
    top_50 = tumor_cells.var[tumor_cells.var['highly_variable']].index[:50].tolist()
    if 'PCSK1N' in tumor_cells.var_names and 'PCSK1N' not in top_50:
        top_50[0] = 'PCSK1N'

    features = tumor_cells[:, top_50].X.toarray() if hasattr(tumor_cells.X, "toarray") else tumor_cells[:, top_50].X
    labels = tumor_cells.obs['Resistance_Status']

    rf_model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
    cv_strategy = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_results = cross_validate(rf_model, features, labels, cv=cv_strategy, scoring=['accuracy', 'f1_macro'], return_estimator=True)

    pd.DataFrame({
        'Metric': ['Mean_Accuracy', 'Mean_F1_Score'],
        'Value': [np.mean(cv_results['test_accuracy']), np.mean(cv_results['test_f1_macro'])]
    }).to_csv(f"{data_dir}Cross_Validation_Metrics.csv", index=False)
    
    rf_model.fit(features, labels)
    pd.DataFrame({
        'Biomarker_Gene': top_50, 
        'RF_Importance_Weight': rf_model.feature_importances_
    }).sort_values(by='RF_Importance_Weight', ascending=False).to_csv(f"{data_dir}Ablated_Feature_Importances.csv", index=False)

    print("3. Computing Unsupervised PAGA Trajectory...")
    tumor_cells.obs['Sample_Origin'] = tumor_cells.obs['Sample_Origin'].astype('category')
    sc.pp.neighbors(tumor_cells, n_neighbors=15, n_pcs=30)
    sc.tl.paga(tumor_cells, groups='Sample_Origin')

    categories = tumor_cells.obs['Sample_Origin'].cat.categories
    connectivity_matrix = tumor_cells.uns['paga']['connectivities'].toarray()
    pd.DataFrame(connectivity_matrix, index=categories, columns=categories).to_csv(f"{data_dir}PAGA_Evolutionary_Connectivity.csv")

    print("--- Phase 6 Complete. Pipeline successfully consolidated! ---")

if __name__ == "__main__":
    main()