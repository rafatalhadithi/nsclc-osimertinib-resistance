import pandas as pd
import requests
import networkx as nx
import os
import io

# PATH CONFIGURATIONS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '../figures/Fig 7'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 08: PPI NETWORK (FIGURE 7) ---")
    
    ml_features_path = os.path.join(DATA_DIR, 'Reviewer_Feature_Importances.csv')
    if not os.path.exists(ml_features_path):
        raise FileNotFoundError(f"Missing {ml_features_path}. Run Script 03 first.")
        
    df = pd.read_csv(ml_features_path)
    genes = df['Gene'].head(50).tolist()
    
    print("Querying STRING Database API...")
    url = "https://version-11-5.string-db.org/api/tsv/network"
    response = requests.post(url, data={"identifiers": "%0d".join(genes), "species": 9606, "required_score": 400})
    
    if response.status_code == 200:
        edges_df = pd.read_csv(io.StringIO(response.text), sep='\t')
        
        G = nx.Graph()
        for _, row in edges_df.iterrows():
            G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
        
        graphml_path = os.path.join(FIG_DIR, 'Figure_7_PPI_Network.graphml')
        nx.write_graphml(G, graphml_path)
        
        print(f"Success! Network with {G.number_of_nodes()} nodes saved to {graphml_path}.")
    print("--- SCRIPT 08 COMPLETE ---")

if __name__ == "__main__":
    main()
