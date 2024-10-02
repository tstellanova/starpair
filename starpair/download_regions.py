import sys
from time import perf_counter

import astropy
from astropy.coordinates import SkyCoord
from astropy.io.fits import VerifyError
from astropy.table import Table
from astroquery.gaia import Gaia
import numpy as np
from astropy import units as u
from astropy.units import Quantity


def collect_one_cone(axis_lon: float, axis_lat: float, cone_radius_deg: float, outfile_prefix=None):
    # axis_coord = SkyCoord(l=180*u.degree, b=0*u.degree, frame='galactic')
    axis_coord = SkyCoord(l=axis_lon * u.degree, b=axis_lat * u.degree, frame='galactic')
    radius = u.Quantity(cone_radius_deg, u.deg)
    # print(f"Searching coord: {axis_coord} \nradius: {cone_radius_deg}")
    int_radius = int(radius.value)
    int_lon = int(axis_lon * 100)
    int_lat = int(axis_lat * 100)

    raw_outfile_name = f"{outfile_prefix}_r{int_radius}_all_{int_lon}_{int_lat}.fits.gz"
    print(f"Searching lon: {axis_lon} lat: {axis_lat} radius: {cone_radius_deg} deg: \n{raw_outfile_name}")

    perf_start_query_region = perf_counter()
    # Perform a cone search on Gaia
    job = Gaia.cone_search_async(coordinate=axis_coord,
                                 radius=radius,
                                 columns=['DESIGNATION', 'ref_epoch', 'ra', 'dec', 'parallax', 'pm', 'pmra', 'pmdec',
                                          'phot_g_mean_mag'],
                                 dump_to_file=True,
                                 output_format='fits',
                                 output_file=raw_outfile_name
                                 )
    star_table = job.get_results()
    print(f"region query took: {(perf_counter() - perf_start_query_region):0.2f} seconds")

    # print(f"star_table  info:")
    # star_table.info()

    parallax = star_table['parallax']  # in milliarcseconds
    # try to filter out stars with unknown parallax
    for i, value in enumerate(parallax):
        # Check if the value is a numpy array, indicating unknown parallax
        if isinstance(value, np.ndarray):
            parallax[i] = sys.float_info.epsilon

    parallax = star_table['parallax']  # in milliarcseconds
    distance_pc = 1000 / parallax  # Convert parallax to distance in parsecs

    # Filter star_table within 30,000 parsecs (~100,000 light years)
    # Fitlter star table to nearby stars (in parsecs)
    max_distance = 400  # 30000  # 30,000 parsecs is a rough boundary for the Milky Way
    # 400 pc ~= 1200 light-years
    max_distance = int(1000/3.26156) # 1000 light-years in parsecs

    # Filter out distant objects
    nearby_stars = star_table[distance_pc <= max_distance]
    # print(f"nearby:  {nearby_stars.info}")
    filtered_outfile_name = f"{outfile_prefix}_r{int_radius}_d{max_distance}_{int_lon}_{int_lat}.fits"
    print(f"nearby stars: {len(nearby_stars)}/{len(star_table)} writing to: \n{filtered_outfile_name}")


    try:
        nearby_stars.write(filtered_outfile_name, format='fits', overwrite=True)
    except VerifyError:
        print(f"Ignoring VerifyError")
        pass


def collect_anti_galactic_cone_series(outfile_prefix=None):
    # subcone_radius_deg = 30
    # cone_axes = ( (180, 0),
    #               (195, 0),
    #               ( 187.5, 13 ),
    #               ( 172.5, 13 ),
    #               ( 165, 0),
    #               ( 172.5, -13 ),
    #               ( 187.5, -13 ),
    #               )

    subcone_radius_deg = 10
    cone_axes = (
        # first layer of 6 cones surrounding center cone:
        (190, 0),
        (185, 8.66),
        (175, 8.66),
        (170, 0),
        (175, -8.66),
        (185, -8.66),
        # second layer of 8 cones surrounding the inner core
        (200, 0),
        (194.14, 14.14),
        (180, 20),
        (165.86, 14.14),
        (160, 0),
        (165.86, -14.14),
        (180, -20),
        (194.14, -14.14),
        # central cone:
        (180, 0)
    )

    for axis_lon, axis_lat in cone_axes:
        collect_one_cone(axis_lon, axis_lat, subcone_radius_deg, outfile_prefix)


# Query a region of the sky
# radius = u.Quantity(30, u.deg)
# galactic_center = SkyCoord(l=0*u.degree, b=0*u.degree, frame='galactic')
# anti_galactic_line = SkyCoord(l=180*u.degree, b=0*u.degree, frame='galactic')
# coord = anti_galactic_line
# coord = galactic_center
Gaia.ROW_LIMIT = 10000  # To return an unlimited number of rows set Gaia.ROW_LIMIT to -1.
# collect_cone_series(anti_galactic_line,outfile_prefix="./data/anti_galactic_line")
# collect_cone_series(gaddlactic_center,outfile_prefix="./data/galactic_line")
max_stars = Gaia.ROW_LIMIT
limit_label_k = int(np.ceil(max_stars/1000))
collect_anti_galactic_cone_series(outfile_prefix=f"./data/antigalactic_L{limit_label_k}k")
