"""
Download Gaia source IDs for Tess TIC IDs 
"""
import pandas as pd
from astroquery.mast import Catalogs
import time

def get_gaia_source_id(tic_id):
    """
    Query MAST for the nearest Gaia DR3 source ID using the TIC ID.
    Returns the Gaia source ID if available, otherwise returns None.
    """
    try:
        catalog_data = Catalogs.query_object(f"TIC {tic_id}", catalog="TIC", version=8)
        if len(catalog_data) > 0:
            # Gaia source ID is stored under the column 'GAIA'
            gaia_id = catalog_data['GAIA'][0]
            return gaia_id
        else:
            return None
    except Exception as e:
        print(f"Error retrieving Gaia ID for TIC {tic_id}: {e}")
        return None

def process_toi_csv(file_path):
    """
    Processes the TESS TOI CSV file, extracts unique TIC IDs, and queries for the nearest Gaia DR3 source IDs.
    """
    # Read the CSV file
    toi_data = pd.read_csv(file_path)

    # Extract unique TIC IDs from the 'TIC' column
    tic_ids = toi_data['TIC'].unique()

    # Create a dictionary to store TIC to Gaia DR3 mappings
    tic_to_gaia = {}

    # Iterate through TIC IDs and retrieve Gaia DR3 source IDs
    for tic_idx in range(len(tic_ids)):
        tic_id = tic_ids[tic_idx]
        gaia_id = get_gaia_source_id(tic_id)
        if gaia_id:
            tic_to_gaia[tic_id] = gaia_id
            print(f"{tic_idx},{tic_id},{gaia_id}")
        else:
            # tic_to_gaia[tic_id] = "No Gaia source ID found"
            print(f"### {tic_idx} {tic_id} -> Invalid")

    time.sleep(0.000001) # Adding a slight delay to avoid hitting API rate limits

    return tic_to_gaia

def save_to_csv(tic_to_gaia, output_file):
    """
    Saves the TIC to Gaia DR3 source ID mapping to a CSV file.
    """
    df = pd.DataFrame(list(tic_to_gaia.items()), columns=['TIC', 'Gaia_DR3_ID'])
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Main script
if __name__ == "__main__":
    input_file = "./tess/2024_10_10_tois.csv"
    output_file = "./tess/tic_to_gaia_mapping.csv"

    tic_to_gaia_mapping = process_toi_csv(input_file)
    save_to_csv(tic_to_gaia_mapping, output_file)