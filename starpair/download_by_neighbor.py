import argparse
from time import perf_counter

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv


#     DESIGNATION,ref_epoch,ra,dec,parallax,l,b,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc,
#     DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep


# Function to query Gaia for stars around a specific source
# def query_gaia_for_source_neighbors(source_id, max_radial_dist_pc: float = 100)

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

# Function to query Gaia for stars within 0.25 degrees and 100 parsecs
# This will batch multiple coordinates into one query
def query_nearby_stars_batch(center_source_id, coordinates, max_dist_pc: float = 100.):
    all_results = []
    for l, b in coordinates:
        # query = f"""
        # SELECT gaia.SOURCE_ID,gaia.l,gaia.b,ABS(1000./gaia.parallax) AS dist_pc,
        #     gaia.distance_gspphot,gaia.parallax,gaia.ra,gaia.dec,gaia.phot_g_mean_mag
        # FROM gaiadr3.gaia_source AS gaia
        # WHERE DISTANCE(POINT('GALACTIC', gaia.l, gaia.b), POINT('GALACTIC', {l}, {b})) < 0.25
        # AND gaia.parallax is not NULL
        # AND gaia.parallax > 0
        # AND 1000.0/gaia.parallax <= {max_dist_pc}
        # """

        query = f"""
        SELECT SOURCE_ID,l,b,  
            distance_gspphot,parallax,ra,dec,phot_g_mean_mag,
            DISTANCE( POINT({l},{b}), POINT(l, b)) AS ang_sep,
            1000./parallax AS dist_pc
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS( POINT({l},{b}), CIRCLE(l, b, 0.25))
        AND SOURCE_ID != {center_source_id}
        AND parallax is not NULL
        AND parallax > 0
        AND 1000.0/parallax <= {max_dist_pc}
        """
        # print(f" neighbor query: \n{query}")
        job = Gaia.launch_job_async(query)
        results = job.get_results()
        all_results.append(results)

    return np.vstack(all_results) if all_results else None


g_habstar_by_id_map = {}
# g_neighbors_by_habstar_id = {}

def main():
    parser = argparse.ArgumentParser(description='Download stars near other stars of interest')
    parser.add_argument('-f', dest='hab_path', nargs="?",
                        default="./data/habstars_2024_10_15.csv",
                        help="csv habitability catalog file",
                        )
    parser.add_argument('-o', dest='outdir', type=str, default='./data/',
                        help='Output directory processed files')

    args = parser.parse_args()
    toi_path = args.hab_path
    out_dir = args.outdir

    total_gaia_objects = 2E9  #overestimate
    Gaia.ROW_LIMIT = int(total_gaia_objects)
    outfile_prefix = f"{out_dir}neighbors_L{total_gaia_objects:0.0e}"
    max_dist_pc = 100
    max_ang_sep = 0.25
    int_ang_sep = 1000 * max_ang_sep
    filtered_outfile_name = f"{outfile_prefix}_a{int_ang_sep}_d{max_dist_pc}.csv"

    print(f"Loading TOIs from: {toi_path} , output to: {filtered_outfile_name}")

    hab_list_row_count = 0
    # import the habitable star system list file
    with open(toi_path, mode='r') as hab_list_file:
        csv_reader = csv.DictReader(hab_list_file)
        while True:
            row_dict = next(csv_reader, None)
            if row_dict is not None:
                hab_list_row_count += 1
                source_id = row_dict['source_id']
                source_id = np.uint64(source_id)
                g_habstar_by_id_map[source_id] = row_dict
            else:
                break
    n_concrete_habitable_stars = hab_list_row_count
    print(f"n_concrete_habitable_stars: {hab_list_row_count}")

    # setup the output file
    field_names = ["orig_id","dest_id","orig_dist","dest_dist","ang_sep","radial_sep","orig_coord","dest_coord"]
    file_ref = open(filtered_outfile_name, 'a')
    csv_writer = csv.writer(file_ref)
    csv_writer.writerow(field_names)
    file_ref.flush()

    evaluated_hstars_count = 0
    total_star_pairs = 0
    for hab_star_source_id in g_habstar_by_id_map.keys():
        habstar_info = g_habstar_by_id_map[hab_star_source_id]
        habstar_l = np.float64(habstar_info['l'])
        habstar_b = np.float64(habstar_info['b'])
        habstar_raw_coord = [(habstar_l, habstar_b)]
        habstar_gspphot_pc = habstar_info['distance_gspphot']
        if habstar_gspphot_pc != '--':
            habstar_radial_dist_pc = np.float64(habstar_gspphot_pc)
        else:
            habstar_radial_dist_pc = np.float64(habstar_info['dist_pc'])
        habstar_coord_str = f"{habstar_l:0.4f} {habstar_b:0.4f} {habstar_radial_dist_pc:0.4f}"

        # Query for all neighbor stars for a single habstar
        neighbors = query_nearby_stars_batch(hab_star_source_id, habstar_raw_coord)
        evaluated_hstars_count += 1
        n_new_neighbors = 0
        if neighbors is not None:
            neighbors = neighbors[0]
            n_new_neighbors = len(neighbors)
            total_star_pairs += n_new_neighbors
            # g_neighbors_by_habstar_id[hab_star_source_id] = neighbors
            for neighbor in neighbors:
                ang_sep = np.float64(neighbor['ang_sep'])
                neighbor_gspphot_pc = neighbor['distance_gspphot']
                if neighbor_gspphot_pc != '--':
                    neighbor_radial_dist_pc = np.float64(neighbor_gspphot_pc)
                else:
                    neighbor_radial_dist_pc = np.float64(neighbor['dist_pc'])

                neighbor_id = np.uint64(neighbor['SOURCE_ID'])
                neighbor_l = np.float64(neighbor['l'])
                neighbor_b = np.float64(neighbor['b'])
                neighbor_coord_str = f"{neighbor_l:0.4f} {neighbor_b:0.4f} {neighbor_radial_dist_pc:0.4f}"

                radial_sep = neighbor_radial_dist_pc - habstar_radial_dist_pc
                if radial_sep > 0:
                    origin_id = neighbor_id
                    dest_id = hab_star_source_id
                    origin_dist = neighbor_radial_dist_pc
                    dest_dist = habstar_radial_dist_pc
                    origin_coord = neighbor_coord_str
                    dest_coord = habstar_coord_str
                else:
                    origin_id = hab_star_source_id
                    dest_id = neighbor_id
                    origin_dist = habstar_radial_dist_pc
                    dest_dist = neighbor_radial_dist_pc
                    origin_coord = habstar_coord_str
                    dest_coord = neighbor_coord_str

                out_row = (origin_id, dest_id, origin_dist, dest_dist, ang_sep, np.abs(radial_sep), origin_coord, dest_coord)
                csv_writer.writerow(out_row)
            file_ref.flush()

        print(f"{evaluated_hstars_count}/{n_concrete_habitable_stars} {hab_star_source_id} has: {n_new_neighbors} neighbors")


    file_ref.flush()
    file_ref.close()
    print(f"Total pairs: {total_star_pairs}")
    print(f"Wrote to: \n{filtered_outfile_name} ")



if __name__ == "__main__":
    main()
