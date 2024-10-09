"""
Search local file Table of stars for star pairs,
where we're looking for stars that are both along an axis "in front" of Earth,
not with Earth between the two stars in the pair.

To implement this we search star tables that are already sorted into cones
either centered on the galactic line (pointed toward the Milky Way galactic center),
or on the anti-galactic line (pointed way from the galactic center).

We first collect all possible star pairs in these tables, then sort them in ascending order by their most distant node.
Star pairs whose midpoint is closer to Earth would likely have more power receivable at Earth.
Secondary sort could be by angular separation between the two stars in the pair,
with the idea that Earth is more likely to receive overshoot signals with shallower glancing angles.


We expect star tables to be stored in FITS files with at least these columns from gaiadr3.gaia_source:

DESIGNATION,ref_epoch,ra,dec,parallax,pm,pmra,pmdec,phot_g_mean_mag,
ABS(1000./parallax) AS dist_pc,
DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep

ORDER BY dist_pc, ang_sep


"""
import csv
import math
import os
import sys
from itertools import combinations
from time import perf_counter

import astropy
from astropy.coordinates import SkyCoord, Angle, Distance
# from astroquery.gaia import Gaia
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.units import Quantity
from astropy.io import fits

import argparse

import functools
from tqdm import tqdm


# import cProfile


# Define a custom cache decorator that caches based only on the first argument
def cache_by_first_arg(func):
    @functools.cache
    def cached_func(first_arg):
        # Return a function that accepts second_arg for actual computation
        def compute_with_second_arg(second_arg):
            return func(first_arg, second_arg)

        return compute_with_second_arg

    def wrapper(first_arg, second_arg):
        # Retrieve the cached function based on the first argument
        compute_with_second_arg = cached_func(first_arg)
        # Call the computation function with the actual second argument
        return compute_with_second_arg(second_arg)

    return wrapper


def star_source_id_from_desig(star_desig: str) -> np.uint64:
    id_val_str = star_desig.replace('Gaia DR3 ', '')
    id_val_long = np.uint64(id_val_str)
    return id_val_long


# g_n_expected_combos = 0
# g_n_evaluated_combos = 0
# g_n_found_pairs = 0

g_starid_star_map = {}


def extract_star_pos(star):
    # pre-calculated distance from observer to star, in parsecs
    dist_pc = float(star['dist_pc'])
    # 	l (galactic longitude) ranges from 0 deg to 360 deg
    galactic_longitude = star['l']
    # 	b (galactic latitude) ranges from -90 deg to +90 deg
    galactic_latitude = star['b']
    # ADQL precalculated: DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep
    axis_sep = star['ang_sep']
    return galactic_longitude, galactic_latitude, dist_pc, axis_sep


g_starid_allcoord_map = {}


@functools.cache
def star_coord_lookup(star_id: np.uint64):
    all_coord = g_starid_allcoord_map.get(star_id)
    if all_coord is None:
        star = g_starid_star_map[star_id]
        lon, lat, rdist, axis_sep = extract_star_pos(star)
        skycoord = SkyCoord(l=lon * u.deg,
                            b=lat * u.deg,
                            distance=rdist,
                            frame='galactic')
        all_coord = lon, lat, rdist, axis_sep, skycoord
        g_starid_allcoord_map[star_id] = all_coord

    return all_coord


def range_filter_ok(lon1: np.float64, lon2: np.float64,
                    lat1: np.float64, lat2: np.float64,
                    max_ang_sep: np.float64) -> bool:
    if np.abs(lon1 - lon2) > max_ang_sep:
        return False
    elif np.abs(lat1 - lat2) > max_ang_sep:
        return False

    return True


g_skycoord_strs_by_id = {}


@cache_by_first_arg
def skycoord_str_lookup(source_id: np.uint64, skycoord: SkyCoord):
    coord_str = g_skycoord_strs_by_id.get(source_id)
    if coord_str is None:
        coord_str = skycoord.to_string('decimal')
        g_skycoord_strs_by_id[source_id] = coord_str

    return coord_str


def evaluate_one_pair_candidate(src_id1: np.uint64, src_id2: np.uint64, max_ang_sep_deg: np.float64):
    lon1, lat1, rdist1, axis_sep1, skycoord1 = star_coord_lookup(src_id1)
    lon2, lat2, rdist2, axis_sep2, skycoord2 = star_coord_lookup(src_id2)

    # Spatial filter: If the lat or lon is very different between the two stars,
    # then we know the angular separation between them is also bigger than what
    # the (expensive) SkyCoord.separation calculation would return below, and we exit early.
    if not range_filter_ok(lon1, lon2, lat1, lat2, max_ang_sep_deg):
        return None

    # print(f"check {src_id1} > {src_id2}")
    ang_sep: float = float(skycoord1.separation(skycoord2).value)
    if ang_sep > max_ang_sep_deg:
        return None

    # Assume that, if the angular separation between the two points is small,
    # the linear separation between them is about the same as the
    # difference between their radial distances from our observation reference point.
    linear_sep = float(abs(rdist1 - rdist2))

    max_node_dist = max(rdist1, rdist2)
    min_node_dist = min(rdist1, rdist2)
    star1_name = str(src_id1)
    star2_name = str(src_id2)

    if rdist1 > rdist2:
        src_name: str = star1_name
        dest_name: str = star2_name
        src_coord_str = skycoord_str_lookup(src_id1, skycoord1)
        dest_coord_str = skycoord_str_lookup(src_id2, skycoord2)
    else:
        src_name: str = star2_name
        dest_name: str = star1_name
        src_coord_str = skycoord_str_lookup(src_id2, skycoord2)
        dest_coord_str = skycoord_str_lookup(src_id1, skycoord1)

    # print(f"found {src_name} > {dest_name}")

    res = (
        max_node_dist, min_node_dist,
        ang_sep, linear_sep,
        src_coord_str, dest_coord_str,
        src_name, dest_name)

    return res


def pre_cache_stars(stars: Table):
    per_start_precache = perf_counter()
    print(f"Pre-caching {len(stars)} stars")
    perf_start_gen_ids = perf_counter()
    for star in stars:
        star_id = star_source_id_from_desig(star['DESIGNATION'])
        g_starid_star_map[star_id] = star
    print(f"Gen IDs >> elapsed: {perf_counter() - perf_start_gen_ids:0.2f} secs")

    perf_start_gen_coords = perf_counter()
    for star_id in g_starid_star_map.keys():
        lon, lat, rdist, axis_sep, skycoord = star_coord_lookup(star_id)
        skycoord_str_lookup(star_id, skycoord)
    print(f"Calc Coords >> elapsed: {perf_counter() - perf_start_gen_coords:0.2f} secs")

    print(f"precache_stars >> elapsed: {perf_counter() - per_start_precache:0.2f} secs")
    pass




def find_close_pairs(stars: Table, max_ang_sep_deg: np.float64 = 0.25):
    # csv_writer=None, file_ref=None):
    """Search all the given stars for very close neighbors,
    where the maximum on-sky angle between the two stars is less than max_angle_deg,
    and the maximum linear (radial) distance is less than max_radial_dist_pc
    """
    num_stars = len(stars)
    n_expected_combos = math.comb(num_stars, 2)

    # A large number of stars results in an explosion of star combinations that need to be checked
    # for appropriate glancing angle. For example, 409,800 stars in a cone from 0 to 90 degrees separation
    # from the galactic centerline results in 83,967,815,100 possible combinations that need to be checked.

    perf_start_combos = perf_counter()

    # profiler = cProfile.Profile()
    # profile_results = []

    pre_cache_stars(stars)
    all_star_ids = g_starid_star_map.keys()

    print(f"Starting map-filter-list on {n_expected_combos} combos...")
    perf_start_mapping = perf_counter()

    close_pairs = list(
        filter(lambda y: y is not None,
               map(lambda x:
                   evaluate_one_pair_candidate(x[0], x[1], max_ang_sep_deg),
                   tqdm(combinations(all_star_ids, 2),total=n_expected_combos)
                   )
               )
    )
    print(f"map-filter-list >> elapsed: {perf_counter() - perf_start_mapping:0.1f} secs")
    print(f"found: {len(close_pairs)} close pairs")

    return close_pairs


def main():
    parser = argparse.ArgumentParser(description='Process stars in FITS file')
    parser.add_argument('src_path', nargs='?',
                        default="./data/galactic_L2e+09_r90_d100_0_0.fits.gz",
                        # default="./data/galactic_L2e+09_r45_d100_0_0.fits.gz",
                        help="Source data file with `.fits` or `.fits.gz` extension",
                        )
    parser.add_argument('-o', dest='outdir', type=str, default='./data',
                        help='Output directory processed files')

    args = parser.parse_args()
    data_file_name = args.src_path
    basename_without_ext = os.path.splitext(os.path.basename(data_file_name))[0]

    out_dir = args.outdir
    print(f"Loading: {data_file_name} , output to: {out_dir}")

    hdulist = fits.open(data_file_name)
    print(f"hdulist: {hdulist}")
    # hdulist[0] is a kind of FITS header/descriptor?
    # hdulist[1] contains the main interesting table
    stars_table = Table.read(hdulist[1])
    print(f"stars_table ({len(stars_table)}): ")
    stars_table.info()
    # print(f"dec units: {stars_table['dec'].unit}")

    # glancing_angle_deg = 0.25
    # max_glancing_angle = Angle(f'{glancing_angle_deg:0.2}d')
    max_glancing_angle: np.float64 = 0.25
    max_radial_dist_pc = 60
    max_glancing_angle_int = int(max_glancing_angle * 1000)

    field_names = ["max_node_dist", "min_node_dist", "ang_sep", "radial_sep", "coord1", "coord2", "origin_name",
                   "dest_name"]

    close_pairs = find_close_pairs(stars_table, max_ang_sep_deg=max_glancing_angle)

    with open(f"./data/{basename_without_ext}_ma{max_glancing_angle_int}_mr{int(max_radial_dist_pc)}_pairs.csv",
              'w') as file_ref:
        csv_writer = csv.writer(file_ref)
        csv_writer.writerow(field_names)
        file_ref.flush()
        # close_pairs = find_close_pairs(stars_table,
        #                                max_ang_sep_deg=max_glancing_angle,
        #                                csv_writer=csv_writer, file_ref=file_ref)
        csv_writer.writerows(close_pairs)
        file_ref.flush()
        file_ref.close()

    # sort pairs in order by ascending distance to most distant node
    print(f"Sorting all pairs... {len(close_pairs)}")
    close_pairs.sort(key=lambda tup: tup[0])
    print(f"Writing {len(close_pairs)} close pairs")
    with open(f"./data/{basename_without_ext}_ma{max_glancing_angle_int}_mr{int(max_radial_dist_pc)}_sorted.csv",
              'w') as f:
        writer = csv.writer(f)
        writer.writerow(field_names)
        writer.writerows(close_pairs)
        f.flush()
        f.close()


if __name__ == "__main__":
    main()
