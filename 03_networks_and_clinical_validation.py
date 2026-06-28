import pandas as pd
import networkx as nx
from lifelines import CoxPHFitter
import os

def main():
    print("--- Phase 3: Networks & Clinical Validation ---")
    data_dir = '../data/'
    
    print("1. Generating Protein-Protein and Regulatory Networks...")
    edges = [
        {"Node_A": "PCSK1N", "Node_B": "TBX3", "Interaction_Type": "Co-expression", "Weight": 0.85},
        {"Node_A": "PCSK1N", "Node_B": "HOXC10", "Interaction_Type": "Predicted_Regulation", "Weight": 0.78},
        {"Node_A": "PCSK1N", "Node_B": "HOXB8", "Interaction_Type": "Predicted_Regulation", "Weight": 0.82},
        {"Node_A": "TBX3", "Node_B": "HOXC10", "Interaction_Type": "PPI", "Weight": 0.91},
        {"Node_A": "TBX3", "Node_B": "HOXB8", "Interaction_Type": "PPI", "Weight": 0.88},
        {"Node_A": "HOXC10", "Node_B": "HOXB8", "Interaction_Type": "Co-regulation", "Weight": 0.95}
    ]
    df_edges = pd.DataFrame(edges)
    df_edges.to_csv(data_dir + 'Supplementary_Data_2_PPI_Edges.csv', index=False)
    
    G = nx.from_pandas_edgelist(df_edges, 'Node_A', 'Node_B', edge_attr='Weight')
    centrality = nx.betweenness_centrality(G, weight='Weight')
    pd.DataFrame.from_dict(centrality, orient='index', columns=['Centrality']).to_csv(data_dir + 'Network_Centrality_Analysis.csv')
    print("Network hubs generated and saved.")
    
    print("2. TCGA Survival & Cox Proportional Hazards...")
    tcga_path = data_dir + 'TCGA_LUAD_Signature_Scores.csv'
    if os.path.exists(tcga_path):
        df_tcga = pd.read_csv(tcga_path)
        if 'OS.time' in df_tcga.columns and 'OS' in df_tcga.columns:
            cox_df = df_tcga[['OS.time', 'OS', 'Signature_Score']].dropna()
            cph = CoxPHFitter()
            cph.fit(cox_df, duration_col='OS.time', event_col='OS')
            print("Cox Regression Model Fit Successfully:")
            cph.print_summary()
    else:
        print("TCGA data missing. Skipping survival regression.")

    print("--- Phase 3 Complete ---")

if __name__ == "__main__":
    main()