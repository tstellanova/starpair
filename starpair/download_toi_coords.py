import argparse
from time import perf_counter

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv



# Function to get the Galactic coordinates (l, b) for a batch of source_ids
def get_star_info_by_gaia_source_id(source_ids):
    source_id_list = ', '.join(map(str, source_ids))
    query = f"""
    SELECT SOURCE_ID,l,b,ABS(1000./parallax) AS dist_pc, distance_gspphot,parallax,ra,dec,phot_g_mean_mag
    FROM gaiadr3.gaia_source
    WHERE SOURCE_ID IN ({source_id_list})
    AND parallax is not NULL
    AND parallax > 0
    """

    # print(f"coord query:\n{query}")
    job = Gaia.launch_job_async(query)
    results = None
    try:
        results =  job.get_results()
    except Exception:
        print(f"Couldn't retrieve results")

    return results


g_habstar_by_id_map = {}
g_neighbors_by_habstar_id = {}

def main():
    parser = argparse.ArgumentParser(description='Download stars near other stars of interest')
    parser.add_argument('-f', dest='hab_path', nargs="?",
                        default="./tess/tic_to_gaia_mapping.csv",
                        help="csv habitability catalog file with rows: external ID, gaia source ID",
                        )
    parser.add_argument('-o', dest='outdir', type=str, default='./data/',
                        help='Output directory processed files')

    args = parser.parse_args()
    toi_path = args.hab_path
    out_dir = args.outdir

    total_gaia_objects = 2E9  #overestimate
    Gaia.ROW_LIMIT = int(total_gaia_objects)
    outfile_prefix = f"{out_dir}habstars"
    filtered_outfile_name = f"{outfile_prefix}_2024_10_15.csv"

    print(f"Loading TOIs from: {toi_path} , output to: {filtered_outfile_name}")

    # Initialize an empty dictionary
    habitable_map = {}
    hab_list_row_count = 0
    # import the habitable star system list file
    with open(toi_path, mode='r') as hab_list_file:
        csv_reader = csv.reader(hab_list_file)
        unused_header = next(csv_reader)
        for hab_list_row in csv_reader:
            hab_list_row_count += 1
            hab_list_uid, gaia_dr3_id, = hab_list_row
            habitable_map[gaia_dr3_id] = hab_list_uid
        hab_list_file.close()
        hab_list_file = None
        hab_list_row = None
    print(f"num raw habitable Gaia sources: {len(habitable_map)}")

    # Batch process the source IDs
    batch_size = 1000  # You can adjust the batch size for optimization
    habstars_id_list = list(habitable_map.keys())
    for i in range(0, len(habstars_id_list), batch_size):
        batch = habstars_id_list[i:i + batch_size]

        # Step 1: Get galactic coordinates for the batch of source IDs
        habstars_batch = get_star_info_by_gaia_source_id(batch)

        if len(habstars_batch) > 0:
            print(f"Got {len(habstars_batch)} habstars for {len(batch)} requested")
            for habstar_info in habstars_batch:
                orig_source_id = np.uint64(habstar_info['SOURCE_ID'])
                g_habstar_by_id_map[orig_source_id] = habstar_info

    n_concrete_habitable_stars = len(g_habstar_by_id_map)
    print(f"n_concrete_habitable_stars: {n_concrete_habitable_stars}")

    # setup the output file
    field_names = ["source_id",
                   "l","b","dist_pc","distance_gspphot","parallax","ra","dec","phot_g_mean_mag",
                   ]
    file_ref = open(filtered_outfile_name, 'w')
    csv_writer = csv.writer(file_ref)
    csv_writer.writerow(field_names)
    file_ref.flush()

    row_count = 0
    for hab_star_source_id in g_habstar_by_id_map.keys():
        habstar_info = g_habstar_by_id_map[hab_star_source_id]
        csv_writer.writerow(habstar_info)
        row_count += 1

    file_ref.flush()
    file_ref.close()
    print(f"Wrote {row_count} rows to: \n{filtered_outfile_name}")


if __name__ == "__main__":
    main()
