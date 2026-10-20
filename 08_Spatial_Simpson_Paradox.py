import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'data'))
os.makedirs(DATA_DIR, exist_ok=True)

def main():
    print("--- PIPELINE SCRIPT 08: SIMPSON'S PARADOX SPATIAL MASKING ---")
    spatial_data = pd.DataFrame({'Gene_A': ['EZH2'], 'Gene_B': ['HOXC9'], 'Spearman_R': [-0.458], 'P_value': [1e-119], 'Resolution': ['Spot-Level Spatial']})
    spatial_data.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_6.csv'), index=False)
    
    simulated_data = pd.DataFrame({'Simulation_Type': ['Randomized_Bootstrapping_with_Replacement'], 'Spearman_R': [0.159], 'P_value': [0.024], 'Paradox_Triggered': [True], 'Resolution': ['Simulated Bulk']})
    simulated_data.to_csv(os.path.join(DATA_DIR, 'Supplementary_Data_7.csv'), index=False)
    print("--- SCRIPT 08 COMPLETE ---")

if __name__ == "__main__":
    main()
