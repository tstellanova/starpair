import argparse
import os
from time import perf_counter

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv

g_habstar_by_id_map = {}


# Function to get the Galactic coordinates (l, b) for a batch of source_ids
def get_star_info_by_gaia_source_id(source_ids):
    source_id_list = ', '.join(map(str, source_ids))
    query = f"""
    SELECT SOURCE_ID,l,b,parallax,(1000.0/parallax) AS dist_pc,ruwe,phot_g_mean_mag
    FROM gaiadr3.gaia_source
    WHERE SOURCE_ID IN ({source_id_list})
    AND parallax is not NULL
    AND parallax > 0
    """


    print(f"coord query:\n{query}")
    job = Gaia.launch_job_async(query)
    results = None
    try:
        results = job.get_results()
    except Exception:
        print(f"Couldn't retrieve results")

    return results


def main():
    parser = argparse.ArgumentParser(description='Download Gaia info for habitable stars')
    parser.add_argument('-f', dest='hab_path', nargs="?",
                        # default="./tess/tess_hab_zone_cat_d10.csv",
                        # default="./tess/tess_hab_zone_cat_d15.csv",
                        # default="./tess/tess_hab_zone_cat_d20.csv",
                        default="./tess/tess_hab_zone_cat_d30.csv",
                        help="csv habitability catalog file with rows containing 'Gaia_ID' fields",
                        )

    args = parser.parse_args()
    input_habcat_path = args.hab_path
    full_path_sans_ext = os.path.splitext(input_habcat_path)[0]
    filtered_outfile_name = f"{full_path_sans_ext}_gaia.csv"

    print(f"Loading hab stars from: {input_habcat_path} , output to: {filtered_outfile_name}")
    Gaia.ROW_LIMIT = int(2E9)  #overestimate

    # Initialize an empty dictionary
    habitable_map = {}
    with open(input_habcat_path, 'r') as input_file:
        reader = csv.DictReader(input_file)
        for row in reader:
            gaia_id = np.uint64(row['Gaia_ID'])
            # note there may be some duplicates, which are ignored
            habitable_map[gaia_id] = gaia_id

    print(f"Nun unique habitable Gaia sources: {len(habitable_map)}")

    # Batch process the source IDs
    batch_size = 1000  # You can adjust the batch size for optimization
    habstars_id_list = list(habitable_map.keys())
    n_habstar_ids = len(habstars_id_list)
    print(f"Num habstar IDs: {n_habstar_ids}")
    n_missing_ids = 0
    for i in range(0, n_habstar_ids, batch_size):
        batch = habstars_id_list[i:i + batch_size]

        # Step 1: Get galactic coordinates for the batch of source IDs
        habstars_batch = get_star_info_by_gaia_source_id(batch)

        if len(habstars_batch) > 0:
            n_missing_ids = n_habstar_ids - len(habstars_batch)
            print(f"Got {len(habstars_batch)} habstars for {len(batch)} requested: {n_missing_ids} missing")
            for habstar_info in habstars_batch:
                orig_source_id = np.uint64(habstar_info['SOURCE_ID'])
                existing_gaia_habstar = g_habstar_by_id_map.get(orig_source_id)
                if existing_gaia_habstar is not None:
                    print(f"!! Duplicate gaia habstar id {orig_source_id} \nwas: {existing_gaia_habstar} \nnew: {habstar_info}")
                g_habstar_by_id_map[orig_source_id] = habstar_info

    n_concrete_habitable_stars = len(g_habstar_by_id_map)
    print(f"n_concrete_habitable_stars: {n_concrete_habitable_stars}")
    if n_missing_ids > 0:
        only_in_thzc_list = habitable_map.keys() - g_habstar_by_id_map.keys()
        print(f"Missing THZC DR2 IDs in DR3: {only_in_thzc_list}")

    # setup the output file
    field_names = ["source_id",
                   "l", "b", "parallax", "dist_pc", "ruwe", "phot_g_mean_mag",
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
