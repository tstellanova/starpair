import argparse
from time import perf_counter

import numpy as np
from astroquery.gaia import Gaia
import csv


# Function to query Gaia for stars within 0.25 degrees and 100 parsecs
# This will batch multiple coordinates into one query
def query_nearby_stars_batch(center_source_id, coordinates, max_dist_pc: float = 100.):
    all_results = []
    for l, b in coordinates:
        query = f"""
        SELECT SOURCE_ID,l,b,  
            distance_gspphot,parallax,ra,dec,phot_g_mean_mag,
            DISTANCE( POINT({l},{b}), POINT(l, b)) AS ang_sep,
            (1000.0/parallax) AS dist_pc
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
g_origstar_by_id_map = {}
g_double_hab_list = []

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
    max_ang_sep = 1.0
    int_ang_sep = int(1000 * max_ang_sep)
    duplicate_output_filename = f"{outfile_prefix}_a{int_ang_sep}_d{max_dist_pc}.csv"
    unique_origins_output_filename = f"{out_dir}uniquepairs_L{total_gaia_objects:0.0e}_a{int_ang_sep}_d{max_dist_pc}.csv"
    double_hab_output_filename = f"{out_dir}doublehab_L{total_gaia_objects:0.0e}_a{int_ang_sep}_d{max_dist_pc}.csv"
    print(f"Loading TOIs from: {toi_path} , output to: {duplicate_output_filename}")

    hab_list_row_count = 0
    # import the habitable star system list file
    with open(toi_path, mode='r') as hab_list_file:
        csv_reader = csv.DictReader(hab_list_file)
        while True:
            row_dict = next(csv_reader, None)
            if row_dict is not None:
                source_id = row_dict['source_id']
                source_id = np.uint64(source_id)
                habstar_radial_dist_pc = np.float64(row_dict['dist_pc'])
                if habstar_radial_dist_pc <= max_dist_pc:
                    g_habstar_by_id_map[source_id] = row_dict
                    hab_list_row_count += 1
                else:
                    print(f"skip habstar {source_id} dist_pc: {habstar_radial_dist_pc}")
            else:
                break
    n_concrete_habitable_stars = hab_list_row_count
    print(f"Num selected habitable stars: {hab_list_row_count}")

    # setup the output file
    field_names = ["orig_id","dest_id","orig_dist","dest_dist","ang_sep","radial_sep","orig_coord","dest_coord"]
    file_ref = open(duplicate_output_filename, 'w')
    csv_writer = csv.writer(file_ref)
    csv_writer.writerow(field_names)
    file_ref.flush()

    perf_start_all_neighbors = perf_counter()
    evaluated_hstars_count = 0
    duplicate_origins_count = 0
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
            for neighbor in neighbors:
                ang_sep = np.float64(neighbor['ang_sep'])
                neighbor_gspphot_pc = neighbor['distance_gspphot']
                neighbor_radial_dist_pc = np.nan
                if neighbor_gspphot_pc != '--':
                    neighbor_radial_dist_pc = np.float64(neighbor_gspphot_pc)
                if np.isnan(neighbor_radial_dist_pc):
                    neighbor_radial_dist_pc = np.float64(neighbor['dist_pc'])
                if np.isnan(neighbor_radial_dist_pc):
                    neighbor_parallax = np.float64(neighbor['parallax'])
                    neighbor_radial_dist_pc = 1000.0/neighbor_parallax

                neighbor_id = np.uint64(neighbor['SOURCE_ID'])
                neighbor_l = np.float64(neighbor['l'])
                neighbor_b = np.float64(neighbor['b'])
                neighbor_coord_str = f"{neighbor_l:0.4f} {neighbor_b:0.4f} {neighbor_radial_dist_pc:0.4f}"

                radial_sep = neighbor_radial_dist_pc - habstar_radial_dist_pc
                if np.isnan(radial_sep):
                    print(f"{neighbor_id} bad sep: {neighbor_radial_dist_pc} - {habstar_radial_dist_pc}")
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
                radial_sep = np.abs(radial_sep)
                out_row = (origin_id, dest_id, origin_dist, dest_dist, ang_sep, radial_sep, origin_coord, dest_coord)
                csv_writer.writerow(out_row)

                double_hab_info = g_habstar_by_id_map.get(neighbor_id)
                if double_hab_info is not None:
                    print(f"double hab! {origin_id} -> {dest_id} radial_sep: {radial_sep} ang_sep: {ang_sep}")
                    g_double_hab_list.append(out_row)

                existing_orig = g_origstar_by_id_map.get(origin_id)
                if existing_orig is None:
                    g_origstar_by_id_map[origin_id] = out_row
                else:
                    duplicate_origins_count += 1
                    existing_radial_sep = existing_orig[5]
                    if radial_sep > 0.5 and radial_sep < existing_radial_sep:
                        g_origstar_by_id_map[origin_id] = out_row
                        print(f"{duplicate_origins_count} dup orig: {origin_id} -> dest: {dest_id} rsep new: {radial_sep} old: {existing_radial_sep}")

            file_ref.flush()

        print(f"{evaluated_hstars_count}/{n_concrete_habitable_stars} {hab_star_source_id} has: {n_new_neighbors} neighbors")

    file_ref.flush()
    file_ref.close()
    file_ref = None
    print(f"Total pairs: {total_star_pairs} duplicates: {duplicate_origins_count}")
    print(f"Neighbor searching >> elapsed: {perf_counter() - perf_start_all_neighbors:0.3f} seconds")
    print(f"Wrote to: \n{duplicate_output_filename} ")

    n_unique_origins = len(g_origstar_by_id_map)
    if n_unique_origins > 0:
        print(f"Total unique origins: {n_unique_origins}")
        # Now write condensed info with unique origin IDs to a separate file
        # field_names = ["orig_id","dest_id","orig_dist","dest_dist","ang_sep","radial_sep","orig_coord","dest_coord"]
        with open(unique_origins_output_filename, mode='w') as file_ref:
            csv_writer = csv.writer(file_ref)
            csv_writer.writerow(field_names)
            file_ref.flush()
            csv_writer.writerows(g_origstar_by_id_map.values())
            file_ref.flush()
            file_ref.close()
            print(f"Wrote unique origins to: \n{unique_origins_output_filename}")

    n_double_hab_pairs = len(g_double_hab_list)
    if n_double_hab_pairs > 0:
        print(f"Have {n_double_hab_pairs} double-hab starpairs!")
        with open(double_hab_output_filename, mode='w') as file_ref:
            # Write double habs into another file
            csv_writer = csv.writer(file_ref)
            csv_writer.writerow(field_names)
            file_ref.flush()
            csv_writer.writerows(g_double_hab_list)
            file_ref.flush()
            file_ref.close()
            print(f"Wrote double hab pairs to: \n{double_hab_output_filename}")

if __name__ == "__main__":
    main()
