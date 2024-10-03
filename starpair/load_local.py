import math
import sys
from itertools import combinations
from time import perf_counter

import astropy
from astropy.coordinates import SkyCoord, Angle
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

    coord = SkyCoord(ra=star['ra'] * u.deg,
                     dec=star['dec'] * u.deg,
                     distance=(1000 / parallax),
                     frame='icrs'
                     )
    return coord


def distance_between_coords(coord1, coord2):
    dist1 = coord1.distance
    dist2 = coord2.distance
    diff = np.abs(dist1 - dist2)
    if diff == 0:
        print(f"zero dist: {coord1} {coord2}")
    return diff


def find_close_pairs(stars: Table, max_angle_deg: Angle = Angle('1.00d')):
    """Search all the given stars for very close neighbors.
    """
    close_pairs = []
    num_items = len(stars)
    n_expected_combos = math.comb(num_items, 2)

    print(f" {type(stars['ra'])}, {type(stars['dec'])}")
    print(stars['ra'].unit)
    print(stars['dec'].unit)
    max_parallax_dist_pc = 306 # 1000 light-years?

    all_star_combos = combinations(stars, 2)
    for star1, star2 in all_star_combos:
        coord1 = star_to_coord(star1)
        coord2 = star_to_coord(star2)
        ang_sep: Angle = coord1.separation(coord2)
        if ang_sep <= max_angle_deg:
            parallax_dist = distance_between_coords(coord1, coord2)
            if 0 < parallax_dist < max_parallax_dist_pc:
                close_pairs.append((ang_sep, parallax_dist, coord1, coord2, star1, star2))
                print(f"{star1['DESIGNATION']} - {star2['DESIGNATION']} ang: {ang_sep}, dist: {parallax_dist} pc")

    print(f"found: {len(close_pairs)} close pairs of {n_expected_combos} total ")
    # print(f"min {min_sep} mean {mean_sep} max {max_sep}")
    return close_pairs


def main():
    parser = argparse.ArgumentParser(description='Process stars in FITS file')
    parser.add_argument('src_path', nargs='?',
                        default="./data/antigalactic_L10k_r10_d306_merged.fits.gz",
                        # default="./data/anti_galactic_line_2_deg.fits.gz",
                        help="Source data file with FITS extension",
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
    print(f"stars_table ({len(stars_table)}: ")
    stars_table.info()
    print(f"dec units: {stars_table['dec'].unit}")
    close_pairs = find_close_pairs(stars_table) # (ang_sep, parallax_dist, coord1, coord2, star1, star2)
    # print(f"close pairs : {close_pairs}")
    close_pairs.sort(key=lambda tup: tup[1])  # sorts in place by ascending distance between pair elemtns
    print(f"total close pairs: {len(close_pairs)}")


if __name__ == "__main__":
    main()
