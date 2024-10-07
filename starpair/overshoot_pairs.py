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
import cProfile

def star_source_id(star):
    return star['DESIGNATION'].replace('Gaia DR3 ', '')


def star_to_coord(star):
    parallax = star['parallax']

    # TODO store an inverse star map of coord to star ?
    if parallax is None:
        print(f"{star['DESIGNATION']} parallax is None")
        parallax = sys.float_info.epsilon
    elif isinstance(parallax, np.ndarray):
        print(f"{star['DESIGNATION']} parallax is bogus array: {parallax}")
        parallax = parallax[0]
        if parallax == 0:
            parallax = sys.float_info.epsilon
    elif isinstance(parallax, (list, tuple)):
        print(f"{star['DESIGNATION']} parallax is list: {parallax}")
        parallax = sys.float_info.epsilon

    # TODO: dist_pc should be precalculated and ready for this
    coord = SkyCoord(l=star['l'] * u.deg,
                     b=star['b'] * u.deg,
                     distance=float(star['dist_pc']),
                     # distance=(1000 / parallax),
                     frame='galactic'
                     )
    return coord


g_star_id_map = {}

#  coord1.to_string('decimal')
def extract_star_pos(star):
    src_id = star_source_id(star)
    # pre-calculated distance from observer to star, in parsecs
    dist_pc = float(star['dist_pc'])
    # 	l (galactic longitude) ranges from 0 deg to 360 deg
    galactic_longitude = star['l']
    # 	b (galactic latitude) ranges from -90 deg to +90 deg
    galactic_latitude = star['b']

    return galactic_longitude, galactic_latitude, dist_pc

@functools.cache
def skycoord_from_raw(lon, lat, dist):
    return SkyCoord(l=lon * u.deg,
             b=lat * u.deg,
             distance=dist,
             frame='galactic')

g_star_skycoord_map = {}
@functools.cache
def calc_star_separation(src_id1, src_id2):
    pos1 = g_star_skycoord_map.get(src_id1)
    if pos1 is None:
        star1 = g_star_id_map[src_id1]
        lon1, lat1, rdist1 = extract_star_pos(star1)
        coord1 = skycoord_from_raw(lon1, lat1, rdist1)
        g_star_skycoord_map[src_id1] = lon1, lat1, rdist1, coord1
    else:
        lon1, lat1, rdist1, coord1 = pos1

    pos2 = g_star_skycoord_map.get(src_id2)
    if pos2 is None:
        star2 = g_star_id_map[src_id2]
        lon2, lat2, rdist2 = extract_star_pos(star2)
        coord2 = skycoord_from_raw(lon2, lat2, rdist2)
        g_star_skycoord_map[src_id2] = lon2, lat2, rdist2, coord2
    else:
        lon2, lat2, rdist2, coord2 = pos2

    linear_sep = float(abs(rdist1 - rdist2))
    ang_sep: float = float(coord1.separation(coord2).value)
    res = (rdist1, rdist2, linear_sep, ang_sep, coord1, coord2)
    # print(f"res: {res}")

    return res

def star_calc_fast(star1, star2):
    src_id1 = star_source_id(star1)
    src_id2 = star_source_id(star2)
    if g_star_id_map.get(src_id1) is None:
        g_star_id_map[src_id1] = star1
    if g_star_id_map.get(src_id2) is None:
        g_star_id_map[src_id2] = star2

    return calc_star_separation(src_id1, src_id2)

def evaluate_one_combo(star1, star2, max_angle_deg:float= 0.25, max_radial_dist_pc: float = 100.):
    # perf_start_eval_one_combo = perf_counter()
    # TODO: We convert stars to coords multiple times, so either optimize the conversion or cache it (or both)
    (rdist1, rdist2, linear_sep, ang_sep, coord1, coord2) = star_calc_fast(star1, star2)
    if ang_sep <= max_angle_deg:
        max_node_dist = max(rdist1, rdist2)
        min_node_dist = min(rdist1, rdist2)
        if max_node_dist <= max_radial_dist_pc:
            star1_name = star_source_id(star1)
            star2_name = star_source_id(star2)
            if rdist1 > rdist2:
                src_name: str = star1_name
                dest_name: str = star2_name
                src_coord = coord1.to_string('decimal')
                dest_coord = coord2.to_string('decimal')
            else:
                src_name: str = star2_name
                dest_name: str = star1_name
                src_coord = coord2.to_string('decimal')
                dest_coord = coord1.to_string('decimal')
            return (
                max_node_dist, min_node_dist,
                ang_sep, linear_sep,
                src_coord, dest_coord,
                src_name, dest_name)
        else:
            return None

            # close_pairs.append(star_pair)
            # n_found_pairs += 1
            # if csv_writer is not None:
            #     csv_writer.writerow(star_pair)
            # print(f"{n_found_pairs} found, {n_evaluated_combos} evaluated, {n_expected_combos} combos "
            #       f">> elapsed : {perf_counter() - perf_start_combos:0.2f} sec")
            # print(f"src: {max_node_dist:0.4f} pc, dst: {min_node_dist:0.4f} pc, "
            #       f"ang_sep: {ang_sep:0.2f}, radial_sep: {linear_sep:0.2f} pc "
            #       f" origin {src_coord} dest {dest_coord} "
            #       f"({src_name} > {dest_name}) ")
            # print(f"One eval >> elapsed: {perf_counter() - perf_start_eval_one_combo:0.4f} sec")
            # print(f"pair: {star_pair}")

def find_close_pairs(stars: Table, max_angle_deg:float= 0.25, max_radial_dist_pc: float = 100.,
                     csv_writer=None):
    """Search all the given stars for very close neighbors,
    where the maximum on-sky angle between the two stars is less than max_angle_deg,
    and the maximum linear (radial) distance is less than max_radial_dist_pc
    """
    close_pairs = []
    num_items = len(stars)
    n_expected_combos = math.comb(num_items, 2)
    n_found_pairs: int = 0
    n_evaluated_combos = 0
    perf_start_combos = perf_counter()
    all_star_combos = combinations(stars, 2)

    profiler = cProfile.Profile()
    profile_results = []

    for star1, star2 in all_star_combos:
        n_evaluated_combos += 1
        profiler.enable()
        found_pair = evaluate_one_combo(star1, star2, max_angle_deg, max_radial_dist_pc)
        profiler.disable()
        profile_results.append(profiler)

        if found_pair is not None:
            n_found_pairs += 1
            close_pairs.append(found_pair)
            if csv_writer is not None:
                csv_writer.writerow(found_pair)
            print(f"found_pair: {found_pair}")

            print(f"{n_found_pairs} found, {n_evaluated_combos} evaluated, {n_expected_combos} combos "
                  f">> elapsed : {perf_counter() - perf_start_combos:0.2f} sec")
            if n_found_pairs % 2 == 0:
                profile_results.print_stats(sort='time')

        #
        # # TODO: We convert stars to coords multiple times, so either optimize the conversion or cache it (or both)
        # (rdist1, rdist2, linear_sep, ang_sep, coord1, coord2) = star_calc_fast(star1, star2)
        # if ang_sep <= max_angle_deg:
        #     max_node_dist = max(rdist1, rdist2)
        #     min_node_dist = min(rdist1, rdist2)
        #     if max_node_dist <= max_radial_dist_pc:
        #         star1_name = star_source_id(star1)
        #         star2_name = star_source_id(star2)
        #         if rdist1 > rdist2:
        #             src_name: str = star1_name
        #             dest_name: str = star2_name
        #             src_coord = coord1.to_string('decimal')
        #             dest_coord = coord2.to_string('decimal')
        #         else:
        #             src_name: str = star2_name
        #             dest_name: str = star1_name
        #             src_coord = coord2.to_string('decimal')
        #             dest_coord = coord1.to_string('decimal')
        #         star_pair = (
        #             max_node_dist, min_node_dist,
        #             ang_sep, linear_sep,
        #             src_coord, dest_coord,
        #             src_name, dest_name)
        #         close_pairs.append(star_pair)
        #         n_found_pairs += 1
        #         if csv_writer is not None:
        #             csv_writer.writerow(star_pair)
        #         print(f"{n_found_pairs} found, {n_evaluated_combos} evaluated, {n_expected_combos} combos "
        #               f">> elapsed : {perf_counter() - perf_start_combos:0.2f} sec")
        #         print(f"src: {max_node_dist:0.4f} pc, dst: {min_node_dist:0.4f} pc, "
        #               f"ang_sep: {ang_sep:0.2f}, radial_sep: {linear_sep:0.2f} pc "
        #               f" origin {src_coord} dest {dest_coord} "
        #               f"({src_name} > {dest_name}) ")
        #         # print(f"One eval >> elapsed: {perf_counter() - perf_start_eval_one_combo:0.4f} sec")
        #         # print(f"pair: {star_pair}")

    print(f"found: {n_found_pairs} ({len(close_pairs)}) close pairs of {n_expected_combos} total ")
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
    max_glancing_angle = 0.25
    max_radial_dist_pc = 60
    max_glancing_angle_int = int(max_glancing_angle * 1000)

    field_names = ["max_node_dist", "min_node_dist", "ang_sep_val", "radial_sep", "coord1", "coord2", "origin_name",
                   "dest_name"]

    with open(f"./data/{basename_without_ext}_ma{max_glancing_angle_int}_mr{int(max_radial_dist_pc)}_pairs.csv",
              'w') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(field_names)
        # csv_writer.writerows(close_pairs)
        find_close_pairs(stars_table,max_angle_deg=max_glancing_angle, max_radial_dist_pc=max_radial_dist_pc, csv_writer=csv_writer)

    # # sort pairs in order by ascending distance to most distant node
    # print(f"Sorting all pairs... {len(close_pairs)}")
    # close_pairs.sort(key=lambda tup: tup[0])
    # print(f"Writing {len(close_pairs)} close pairs")
    # with open(f"./data/{basename_without_ext}_ma{max_glancing_angle_int}_mr{int(max_radial_dist_pc)}_sorted.csv",
    #           'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(field_names)
    #     writer.writerows(close_pairs)


if __name__ == "__main__":
    main()
