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


def star_to_coord(star):
    parallax = star['parallax']
    star_name = star['DESIGNATION']

    if parallax is None:
        print(f"{star_name} parallax is None")
        parallax = sys.float_info.epsilon
    elif isinstance(parallax, np.ndarray):
        # print(f"{star_name} parallax is bogus array: {parallax}")
        parallax = parallax[0]
        if parallax == 0:
            parallax = sys.float_info.epsilon
    elif isinstance(parallax, (list, tuple)):
        print(f"parallax is list: {parallax}")
        parallax = sys.float_info.epsilon

    coord = SkyCoord(l=star['l'] * u.deg,
                 b=star['b'] * u.deg,
                 distance=(1000 / parallax),
                 frame='galactic'
                 )
    # coord = SkyCoord(ra=star['ra'] * u.deg,
    #                  dec=star['dec'] * u.deg,
    #                  distance=(1000 / parallax),
    #                  frame='icrs'
    #                  )
    return coord


def distance_between_coords(coord1: SkyCoord, coord2: SkyCoord) -> float:
    assert coord1.frame.name == 'galactic'
    distance:Distance = coord1.separation_3d(coord2)
    return float(distance.value)


def find_close_pairs(stars: Table, max_angle_deg: Angle = Angle('0.25d'), max_radial_dist_pc:float = 100.):
    """Search all the given stars for very close neighbors.
    """
    close_pairs = []
    num_items = len(stars)
    n_expected_combos = math.comb(num_items, 2)
    n_found_pairs:int = 0

    all_star_combos = combinations(stars, 2)
    for star1, star2 in all_star_combos:
        coord1 = star_to_coord(star1)
        coord2 = star_to_coord(star2)
        ang_sep: Angle = coord1.separation(coord2)
        if ang_sep <= max_angle_deg:
            star1_dist = float(star1['dist_pc'])
            star2_dist = float(star2['dist_pc'])
            max_node_dist = max(star1_dist, star2_dist)
            min_node_dist = min(star1_dist, star2_dist)
            if max_node_dist <= max_radial_dist_pc:
                radial_sep = distance_between_coords(coord1, coord2)
                star1_name = star1['DESIGNATION'].replace('Gaia DR3 ','')
                star2_name = star2['DESIGNATION'].replace('Gaia DR3 ','')
                ang_sep_val = float(ang_sep.value)
                if star1_dist > star2_dist:
                    origin_name:str = star1_name
                    dest_name:str = star2_name
                    origin_coord = coord1.to_string('decimal')
                    dest_coord = coord2.to_string('decimal')
                else:
                    origin_name:str = star2_name
                    dest_name:str = star1_name
                    origin_coord = coord2.to_string('decimal')
                    dest_coord = coord1.to_string('decimal')

                star_pair = (max_node_dist, min_node_dist, ang_sep_val, radial_sep, origin_coord, dest_coord, origin_name, dest_name)
                close_pairs.append((max_node_dist, min_node_dist, ang_sep_val, radial_sep, origin_coord, dest_coord, origin_name, dest_name))
                n_found_pairs += 1
                print(f"{n_found_pairs} origin: {max_node_dist:0.4f} pc, dest: {min_node_dist:0.4f} pc, "
                      f"ang_sep: {ang_sep:0.2f}, radial_sep: {radial_sep:0.2f} pc "
                      f" origin {origin_coord} dest {dest_coord} "
                      f"({origin_name} > {dest_name}) ")
                # print(f"pair: {star_pair}")

    print(f"found: {n_found_pairs} ({len(close_pairs)}) close pairs of {n_expected_combos} total ")
    return close_pairs


def main():
    parser = argparse.ArgumentParser(description='Process stars in FITS file')
    parser.add_argument('src_path', nargs='?',
                        default="./data/galactic_L2e+09_r45_d100_0_0.fits.gz",
                        # default="./data/antigalactic_L10k_r10_d306_merged.fits.gz",
                        # default="./data/anti_galactic_line_2_deg.fits.gz",
                        help="Source data file with `.fits` or `.fits.gz` extension",
                        )
    parser.add_argument('-o', dest='outdir', type=str, default='./data',
                        help='Output directory processed files')

    args = parser.parse_args()
    data_file_name = args.src_path
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

    glancing_angle_deg = 0.25
    max_glancing_angle = Angle(f'{glancing_angle_deg:0.2}d')
    close_pairs = find_close_pairs(stars_table, max_angle_deg=max_glancing_angle, max_radial_dist_pc=60)
    field_names = ["max_node_dist", "min_node_dist", "ang_sep_val", "radial_sep", "coord1", "coord2", "origin_name", "dest_name"]
    max_glancing_angle_int = int(glancing_angle_deg * 1000)

    # sort pairs in order by ascending distance to most distant node
    print(f"Sorting all pairs... {len(close_pairs)}")
    close_pairs.sort(key=lambda tup: tup[0])
    print(f"Writing {len(close_pairs)} close pairs")
    with open(f"./data/galactic_L2e+09_r45_d100_0_0_a{max_glancing_angle_int}_pairs.csv", 'w') as f:
        writer = csv.writer(f)
        writer.writerow(field_names)
        writer.writerows(close_pairs)



if __name__ == "__main__":
    main()
