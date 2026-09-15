import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# 1. PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures/Fig 5'))
os.makedirs(FIG_DIR, exist_ok=True)

def robust_patient_cv(X, y, groups, n_splits=5, n_features=50):
    """Executes Nested Patient-Level GroupKFold CV to prevent data leakage."""
    n_splits = min(n_splits, len(np.unique(groups)))
    gkf = GroupKFold(n_splits=n_splits)
    
    cv_metrics = {'accuracy': [], 'f1_macro': [], 'precision_macro': [], 'recall_macro': []}
    fold_importances = [] 
    
    for train_idx, test_idx in gkf.split(X, y, groups=groups):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        # Nested Feature Selection (Only looks at training data to prevent global leakage)
        feature_selector = SelectKBest(score_func=f_classif, k=n_features)
        X_train_selected = feature_selector.fit_transform(X_train, y_train)
        X_test_selected = feature_selector.transform(X_test)
        
        selected_mask = feature_selector.get_support()
        selected_genes = X.columns[selected_mask]
        
        # Model Training
        rf_model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
        rf_model.fit(X_train_selected, y_train)
        predictions = rf_model.predict(X_test_selected)
        
        # Save feature importances for this fold
        importances = pd.Series(rf_model.feature_importances_, index=selected_genes)
        fold_importances.append(importances)
        
        # Scoring
        cv_metrics['accuracy'].append(accuracy_score(y_test, predictions))
        cv_metrics['f1_macro'].append(f1_score(y_test, predictions, average='macro'))
        cv_metrics['precision_macro'].append(precision_score(y_test, predictions, average='macro', zero_division=0))
        cv_metrics['recall_macro'].append(recall_score(y_test, predictions, average='macro', zero_division=0))

    mean_metrics = pd.DataFrame(cv_metrics).mean().to_frame(name='Mean_Value')
    importance_df = pd.concat(fold_importances, axis=1).fillna(0)
    mean_importances = importance_df.mean(axis=1).sort_values(ascending=False).to_frame(name='Mean_Importance')
    
    return mean_metrics, mean_importances

def main():
    print("--- PIPELINE SCRIPT 03: ML CLASSIFIER & EXTERNAL VALIDATION ---")
    
    # ---------------------------------------------------------
    # PART A: DISCOVERY COHORT (GSE131907) CV & FEATURE RANKING
    # ---------------------------------------------------------
    disc_path = os.path.join(DATA_DIR, 'GSE131907_FINAL_ANNOTATED.h5ad')
    if not os.path.exists(disc_path):
        raise FileNotFoundError(f"Missing discovery dataset: {disc_path}")
        
    print("1. Loading Discovery Cohort...")
    adata_disc = sc.read_h5ad(disc_path)
    
    # Isolate Epithelial cells and valid origins
    tumor_cells = adata_disc[adata_disc.obs['cell_type'] == 'Epithelial cells'].copy()
    valid_origins = ['tLung', 'tL/B', 'PE', 'mLN']
    tumor_cells = tumor_cells[tumor_cells.obs['Sample_Origin'].isin(valid_origins)].copy()
    
    def assign_status(origin):
        if origin in ['tLung', 'tL/B']: return 0 # Primary_Sensitive
        elif origin in ['PE', 'mLN']: return 1   # Metastatic_Resistant
        return -1
    
    tumor_cells.obs['Target'] = tumor_cells.obs['Sample_Origin'].apply(assign_status)
    tumor_cells = tumor_cells[tumor_cells.obs['Target'] != -1].copy()
    
    print("2. Preparing matrices and running GroupKFold CV...")
    X_disc = pd.DataFrame(tumor_cells.X.toarray() if hasattr(tumor_cells.X, "toarray") else tumor_cells.X, 
                          index=tumor_cells.obs.index, columns=tumor_cells.var_names)
    y_disc = tumor_cells.obs['Target']
    groups_disc = tumor_cells.obs['Sample'] # Patient ID column
    
    metrics, importances = robust_patient_cv(X_disc, y_disc, groups_disc, n_splits=5, n_features=50)
    
    print("\n[Discovery CV Metrics]")
    print(metrics)
    metrics.to_csv(os.path.join(DATA_DIR, 'Reviewer_Patient_Level_CV_Metrics.csv'))
    
    importances.index.name = 'Gene'
    importances.reset_index().to_csv(os.path.join(DATA_DIR, 'Reviewer_Feature_Importances.csv'), index=False)
    
    # ---------------------------------------------------------
    # PART B: GENERATING FIGURE 5 (FEATURE IMPORTANCE BARPLOT)
    # ---------------------------------------------------------
    print("\n3. Generating Figure 5...")
    top_ml = importances.head(10).reset_index()
    top_ml.columns = ['Gene', 'Importance_Score']
    
    plt.figure(figsize=(7, 6))
    sns.set_theme(style="ticks", context="paper")
    sns.barplot(data=top_ml, x='Importance_Score', y='Gene', hue='Gene', palette='mako', edgecolor='black', legend=False)
    sns.despine()
    plt.title('Random Forest Feature Importance (Patient-Level CV)', weight='bold')
    plt.xlabel('Importance Score (Mean Decrease Impurity)', weight='bold')
    plt.ylabel('')
    plt.tight_layout()
    
    fig5_pdf = os.path.join(FIG_DIR, 'Figure_5_ML_Features.pdf')
    plt.savefig(fig5_pdf, format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    print(f"   -> Saved: {fig5_pdf}")

    # ---------------------------------------------------------
    # PART C: EXTERNAL VALIDATION (MAYNARD COHORT)
    # ---------------------------------------------------------
    val_path = os.path.join(DATA_DIR, 'maynard_annotated.h5ad')
    if os.path.exists(val_path):
        print("\n4. Loading External Maynard Cohort for Un-leaked Validation...")
        adata_val = sc.read_h5ad(val_path)
        
        # Filter for TN (0) vs PD (1)
        val_sub = adata_val[adata_val.obs['Clinical_Stage'].isin(['TN', 'PD'])].copy()
        val_sub.obs['Target'] = val_sub.obs['Clinical_Stage'].map({'TN': 0, 'PD': 1})
        
        # Intersect top 50 features
        top_50_genes = importances.head(50).index.tolist()
        shared_genes = [g for g in top_50_genes if g in val_sub.var_names]
        print(f"   -> Found {len(shared_genes)}/50 signature genes in Maynard cohort.")
        
        if len(shared_genes) > 0:
            # Train final model on FULL discovery set using shared genes
            X_train_final = X_disc[shared_genes]
            y_train_final = y_disc
            
            # Scale data independently to mitigate baseline batch shifts
            scaler_train = StandardScaler().fit(X_train_final)
            X_train_scaled = scaler_train.transform(X_train_final)
            
            final_rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
            final_rf.fit(X_train_scaled, y_train_final)
            
            # Prepare Maynard Test Data
            X_test_final = pd.DataFrame(val_sub[:, shared_genes].X.toarray() if hasattr(val_sub.X, "toarray") else val_sub[:, shared_genes].X, columns=shared_genes)
            y_test_final = val_sub.obs['Target'].values
            
            scaler_test = StandardScaler().fit(X_test_final)
            X_test_scaled = scaler_test.transform(X_test_final)
            
            # Predict and Score
            test_preds = final_rf.predict(X_test_scaled)
            test_probs = final_rf.predict_proba(X_test_scaled)[:, 1]
            
            val_acc = accuracy_score(y_test_final, test_preds)
            val_auc = roc_auc_score(y_test_final, test_probs)
            
            print("\n[External Validation Metrics (Maynard TN vs PD)]")
            print(f"   -> Accuracy: {val_acc:.4f}")
            print(f"   -> ROC-AUC:  {val_auc:.4f}")
            
            # Save results
            ext_metrics = pd.DataFrame({'Metric': ['Accuracy', 'ROC-AUC'], 'Value': [val_acc, val_auc]})
            ext_metrics.to_csv(os.path.join(DATA_DIR, 'Maynard_External_Validation_Metrics.csv'), index=False)
            
    print("\n--- SCRIPT 03 COMPLETE ---")

if __name__ == "__main__":
    main()
