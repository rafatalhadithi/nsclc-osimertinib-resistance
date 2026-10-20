import pandas as pd
import requests
import networkx as nx
import os
import io
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'figures', 'Fig 8'))
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 09: PPI NETWORK ---")
    ml_features_path = os.path.join(DATA_DIR, 'Supplementary_Data_5.csv')
    if not os.path.exists(ml_features_path): sys.exit(1)
        
    df = pd.read_csv(ml_features_path)
    genes = df[df.columns[0]].head(50).tolist()
    
    response = requests.post("https://version-11-5.string-db.org/api/tsv/network", data={"identifiers": "%0d".join(genes), "species": 9606, "required_score": 400})
    if response.status_code == 200:
        edges_df = pd.read_csv(io.StringIO(response.text), sep='\t')
        G = nx.Graph()
        for _, row in edges_df.iterrows(): G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
        nx.write_graphml(G, os.path.join(FIG_DIR, 'Figure_8_PPI_Network.graphml'))
    print("--- SCRIPT 09 COMPLETE ---")

if __name__ == "__main__":
    main()
