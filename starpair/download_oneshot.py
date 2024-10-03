import sys
from astroquery.utils.tap.model.job import Job
from time import perf_counter

import astropy
from astropy.coordinates import SkyCoord
from astropy.io.fits import VerifyError
from astropy.table import Table
from astroquery.gaia import Gaia
import numpy as np
from astropy import units as u
from astropy.units import Quantity


def collect_one_cone(axis_lon: float, axis_lat: float, cone_radius_deg: float,  max_dist_pc:float=306, outfile_prefix=None):
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

def oneshot_adql_query(axis_lon: float, axis_lat: float, cone_radius_deg: float, max_dist_pc: float, outfile_prefix=None):

    int_radius = int(cone_radius_deg)
    int_lon = int(axis_lon * 100)
    int_lat = int(axis_lat * 100)

    row_limit =  Gaia.ROW_LIMIT
    filtered_outfile_name = f"{outfile_prefix}_r{int_radius}_d{max_dist_pc}_{int_lon}_{int_lat}.fits.gz"

    #     TOP {row_limit}
    query1 = f"""SELECT
    DESIGNATION,ref_epoch,ra,dec,parallax,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc, DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep
    FROM gaiadr3.gaia_source
    WHERE parallax is not NULL
    AND parallax > 0
    AND (1000./parallax) < {max_dist_pc:0.2f}
    AND DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) < {cone_radius_deg:0.2f}
    AND random_index < {row_limit}
    """

    print(f"Sending query1:")
    print(query1)

    perf_start_query1 = perf_counter()
    job1:Job = Gaia.launch_job_async(query1)
    # print(f"job info: {job1}")
    print(f"query1 took: {(perf_counter() - perf_start_query1):0.2f} seconds")

    job1_id = job1.jobid
    query2 = f"""SELECT *
    FROM job_upload."job{job1_id}"
    ORDER BY dist_pc, ang_sep
    """

    print(f"Sending query2:")
    print(query2)
    perf_start_query2 = perf_counter()
    job2:Job = Gaia.launch_job_async(query2,
                                     dump_to_file=True,
                                     output_file=filtered_outfile_name,
                                     output_format=  "fits"
                                     )
    # print(f"job2 info: {job2}")

    star_table = job2.get_results()
    print(f"query2 took: {(perf_counter() - perf_start_query2):0.2f} seconds")
    print(f"nearby stars {len(star_table)} written to: \n{filtered_outfile_name}")



def collect_antigalactic_cone(subcone_radius_deg=10, outfile_prefix="./data/antigalactic"):
    limit_stars_float = float(Gaia.ROW_LIMIT)
    outfile_prefix=f"./data/antigalactic_L{limit_stars_float:0.0e}"
    # subcone_radius_deg = 45 # we'll collect a cone with 90 degree spread across antigalactic hemisphere
    axis_lon, axis_lat = (180, 0)
    oneshot_adql_query(axis_lon, axis_lat, subcone_radius_deg, max_dist_pc=100, outfile_prefix=outfile_prefix)

def collect_galactic_cone(subcone_radius_deg=10, outfile_prefix="./data/galactic"):
    limit_stars_float = float(Gaia.ROW_LIMIT)
    outfile_prefix=f"./data/galactic_L{limit_stars_float:0.0e}"
    # subcone_radius_deg = 45 # we'll collect a cone with 90 degree spread across galactic hemisphere
    axis_lon, axis_lat = (0, 0)
    oneshot_adql_query(axis_lon, axis_lat, subcone_radius_deg, max_dist_pc=100, outfile_prefix=outfile_prefix)

def main():
    # To return an unlimited number of rows set Gaia.ROW_LIMIT to -1.
    total_gaia_objects = 2E9 #overestimate
    Gaia.ROW_LIMIT = int(total_gaia_objects)
    # To return an unlimited number of rows set Gaia.ROW_LIMIT to -1.
    max_stars = Gaia.ROW_LIMIT
    # limit_stars_float = max_stars/1.0
    perf_start_total_query = perf_counter()
    collect_antigalactic_cone(subcone_radius_deg=30.0)
    collect_galactic_cone(subcone_radius_deg=30.0)
    # collect_antigalactic_cone(outfile_prefix=f"./data/antigalactic_L{limit_stars_float:0.0e}")
    # collect_galactic_cone(outfile_prefix=f"./data/galactic_L{limit_stars_float:0.0e}")
    print(f"total {max_stars} query >>> elapsed: {(perf_counter() - perf_start_total_query):0.2f} seconds")

# note that 20000000 ("_L2e+07_") takes about 352-424 seconds on query1
# note that 200000000 ("_L2e+08_") takes about 1192 seconds on query1
# ./data/galactic_L2e+08_r45_d100_0_0.fits.gz


# nearby stars 26742 written to:
# ./data/galactic_L2e+08_r45_d100_0_0.fits.gz
# total 200000000 query >>> elapsed: 963.44 seconds


# 8x as many near-galactic lines as antigalactic:
# nearby stars 23043 written to:
# ./data/antigalactic_L2e+09_r30_d100_18000_0.fits.gz
# nearby stars 185819 written to:
# ./data/galactic_L2e+09_r30_d100_0_0.fits.gz
# total 2x2000000000 query >>> elapsed: 1259.19 seconds

if __name__ == "__main__":
    main()
