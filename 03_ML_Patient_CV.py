import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 5'))
os.makedirs(FIG_DIR, exist_ok=True)

def robust_patient_cv(X, y, groups, n_splits=5, n_features=50):
    unique_groups = np.unique(groups.dropna())
    n_splits = max(2, min(n_splits, len(unique_groups)))
    gkf = GroupKFold(n_splits=n_splits)
    cv_metrics = {'accuracy': [], 'f1_macro': [], 'precision_macro': [], 'recall_macro': []}
    fold_importances = [] 
    
    for train_idx, test_idx in gkf.split(X, y, groups=groups):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        fs = SelectKBest(score_func=f_classif, k=n_features)
        X_train_sel = fs.fit_transform(X_train, y_train)
        X_test_sel = fs.transform(X_test)
        selected_genes = X.columns[fs.get_support()]
        
        rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
        rf.fit(X_train_sel, y_train)
        preds = rf.predict(X_test_sel)
        
        fold_importances.append(pd.Series(rf.feature_importances_, index=selected_genes))
        cv_metrics['accuracy'].append(accuracy_score(y_test, preds))
        cv_metrics['f1_macro'].append(f1_score(y_test, preds, average='macro'))
        cv_metrics['precision_macro'].append(precision_score(y_test, preds, average='macro', zero_division=0))
        cv_metrics['recall_macro'].append(recall_score(y_test, preds, average='macro', zero_division=0))

    mean_metrics = pd.DataFrame(cv_metrics).mean().to_frame(name='Mean_Value')
    mean_importances = pd.concat(fold_importances, axis=1).fillna(0).mean(axis=1).sort_values(ascending=False).to_frame(name='Mean_Importance')
    return mean_metrics, mean_importances

def main():
    print("--- PIPELINE SCRIPT 03: ML CLASSIFIER & EXTERNAL VALIDATION ---")
    disc_path = os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad')
    if not os.path.exists(disc_path): sys.exit(1)
        
    adata_disc = sc.read_h5ad(disc_path)
    cell_type_col = next((c for c in adata_disc.obs.columns if 'cell' in c.lower()), 'cell_type')
    origin_col = next((c for c in adata_disc.obs.columns if 'origin' in c.lower() or 'site' in c.lower()), 'Sample_Origin')
    sample_col = next((c for c in adata_disc.obs.columns if 'patient' in c.lower() or 'sample' in c.lower()), 'Sample')

    tumor_cells = adata_disc[adata_disc.obs[cell_type_col].astype(str).str.contains('Epithelial|Tumor', case=False, na=False)].copy()
    tumor_cells = tumor_cells[tumor_cells.obs[origin_col].isin(['tLung', 'tL/B', 'PE', 'mLN'])].copy()
    
    tumor_cells.obs['Target'] = tumor_cells.obs[origin_col].map({'tLung':0, 'tL/B':0, 'PE':1, 'mLN':1})
    tumor_cells = tumor_cells[tumor_cells.obs['Target'].notna()].copy()
    
    X_disc = pd.DataFrame(tumor_cells.X.toarray() if hasattr(tumor_cells.X, "toarray") else tumor_cells.X, index=tumor_cells.obs.index, columns=tumor_cells.var_names)
    metrics, importances = robust_patient_cv(X_disc, tumor_cells.obs['Target'], tumor_cells.obs[sample_col])
    
    metrics.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_4.csv'))
    importances.index.name = 'Gene'
    importances.reset_index().to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_5.csv'), index=False)
    
    top_ml = importances.head(10).reset_index()
    top_ml.columns = ['Gene', 'Importance_Score']
    plt.figure(figsize=(7, 6))
    sns.barplot(data=top_ml, x='Importance_Score', y='Gene', hue='Gene', palette='mako', legend=False)
    sns.despine()
    plt.title('Random Forest Feature Importance (Patient-Level CV)', weight='bold')
    plt.xlabel('Importance Score (Mean Decrease Impurity)', weight='bold'); plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_5_ML_Features.pdf'), format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    
    print("--- SCRIPT 03 COMPLETE ---")

if __name__ == "__main__":
    main()
