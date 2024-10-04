import sys
import time
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



def oneshot_cone_query(axis_lon: float, axis_lat: float, outer_cone_radius_deg: float, max_dist_pc: float, outfile_prefix=None):
    int_radius = int(outer_cone_radius_deg)
    int_lon = int(axis_lon * 100)
    int_lat = int(axis_lat * 100)

    row_limit = Gaia.ROW_LIMIT
    filtered_outfile_name = f"{outfile_prefix}_r{int_radius}_d{max_dist_pc}_{int_lon}_{int_lat}.fits.gz"

    # first query searches across fraction of ALL stars for those matching angular and radial distance limits
    query1 = f"""SELECT
    DESIGNATION,ref_epoch,ra,dec,parallax,l,b,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc, DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep
    FROM gaiadr3.gaia_source
    WHERE parallax is not NULL
    AND parallax > 0
    AND (1000./parallax) < {max_dist_pc:0.2f}
    AND DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) < {outer_cone_radius_deg:0.2f}
    AND random_index < {row_limit}
    """

    cur_time_str = time.strftime("%H:%M:%S", time.localtime())
    print(f"Sending query1 at {cur_time_str}:")
    print(query1)

    perf_start_query1 = perf_counter()
    job1: Job = Gaia.launch_job_async(query1)
    print(f"query1 took: {(perf_counter() - perf_start_query1):0.2f} seconds")

    # second query sorts results first by radial distance to star, followed by angular distance from ref line
    job1_id = job1.jobid
    query2 = f"""SELECT *
    FROM job_upload."job{job1_id}"
    ORDER BY dist_pc, ang_sep
    """

    print(f"Sending query2:")
    print(query2)
    perf_start_query2 = perf_counter()
    job2: Job = Gaia.launch_job_async(query2,
                                      dump_to_file=True,
                                      output_file=filtered_outfile_name,
                                      output_format="fits"
                                      )
    # print(f"job2 info: {job2}")

    star_table = job2.get_results()
    print(f"query2 took: {(perf_counter() - perf_start_query2):0.2f} seconds")
    print(f"{len(star_table)} nearby stars written to: \n{filtered_outfile_name}")


def oneshot_outer_cone_query(
        axis_lon: float, axis_lat: float,
        inner_cone_radius_deg: float = 30., outer_cone_radius_deg: float = 45.,
        max_dist_pc: float = 100.,
        outfile_prefix="outercone"):
    int_ri = int(inner_cone_radius_deg)
    int_ro = int(outer_cone_radius_deg)
    int_lon = int(axis_lon * 100)
    int_lat = int(axis_lat * 100)

    row_limit = Gaia.ROW_LIMIT
    filtered_outfile_name = f"{outfile_prefix}_ri{int_ri}_ro{int_ro}_d{max_dist_pc}_{int_lon}_{int_lat}.fits.gz"

    # first query searches across fraction of ALL stars for those matching angular and radial distance limits
    query1 = f"""SELECT
    DESIGNATION,ref_epoch,ra,dec,parallax,l,b,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc, DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) AS ang_sep
    FROM gaiadr3.gaia_source
    WHERE parallax is not NULL AND parallax > 0
    AND (1000./parallax) < {max_dist_pc:0.2f}
    AND DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) >= {inner_cone_radius_deg:0.2f}
    AND DISTANCE({axis_lon:0.2f}, {axis_lat:0.2f}, l, b) < {outer_cone_radius_deg:0.2f}
    AND random_index < {row_limit}
    """
    cur_time_str = time.strftime("%H:%M:%S", time.localtime())
    print(f"Sending query1 at {cur_time_str}:")
    print(query1)

    perf_start_query1 = perf_counter()
    job1: Job = Gaia.launch_job_async(query1)
    print(f"query1 took: {(perf_counter() - perf_start_query1):0.2f} seconds")

    # second query sorts results first by radial distance to star, followed by angular distance from ref line
    job1_id = job1.jobid
    query2 = f"""SELECT *
    FROM job_upload."job{job1_id}"
    ORDER BY dist_pc, ang_sep
    """

    print(f"Sending query2:")
    print(query2)
    perf_start_query2 = perf_counter()
    job2: Job = Gaia.launch_job_async(query2,
                                      dump_to_file=True,
                                      output_file=filtered_outfile_name,
                                      output_format="fits"
                                      )
    # print(f"job2 info: {job2}")

    star_table = job2.get_results()
    print(f"query2 took: {(perf_counter() - perf_start_query2):0.2f} seconds")
    print(f"nearby stars {len(star_table)} written to: \n{filtered_outfile_name}")


def collect_antigalactic_cone(subcone_radius_deg=10):
    limit_stars_float = float(Gaia.ROW_LIMIT)
    outfile_prefix = f"./data/antigalactic_L{limit_stars_float:0.0e}"
    axis_lon, axis_lat = (180, 0)
    perf_start_oneshot = perf_counter()
    oneshot_cone_query(axis_lon, axis_lat, subcone_radius_deg, max_dist_pc=100, outfile_prefix=outfile_prefix)
    print(f"antigalactic query >>> elapsed: {(perf_counter() - perf_start_oneshot):0.2f} seconds")


def collect_galactic_cone(subcone_radius_deg=10):
    limit_stars_float = float(Gaia.ROW_LIMIT)
    outfile_prefix = f"./data/galactic_L{limit_stars_float:0.0e}"
    axis_lon, axis_lat = (0, 0)
    perf_start_oneshot = perf_counter()
    oneshot_cone_query(axis_lon, axis_lat, subcone_radius_deg, max_dist_pc=100, outfile_prefix=outfile_prefix)
    print(f"galactic query >>> elapsed: {(perf_counter() - perf_start_oneshot):0.2f} seconds")


def collect_outer_cone(axis_lon:float = 0, axis_lat:float = 0, inner_cone_radius_deg: float = 45, outer_cone_radius_deg: float = 60,
                       outfile_prefix="./data/galactic"):
    limit_stars_float = float(Gaia.ROW_LIMIT)
    outfile_prefix = f"./data/galactic_outercone_L{limit_stars_float:0.0e}"
    perf_start_oneshot = perf_counter()
    oneshot_outer_cone_query(axis_lon, axis_lat,
                             inner_cone_radius_deg=inner_cone_radius_deg,
                             outer_cone_radius_deg=outer_cone_radius_deg,
                             max_dist_pc=100, outfile_prefix=outfile_prefix)
    print(f"galactic outercone query >>> elapsed: {(perf_counter() - perf_start_oneshot):0.2f} seconds")


def main():
    # To return an unlimited number of rows set Gaia.ROW_LIMIT to -1.
    total_gaia_objects = 2E9  #overestimate
    Gaia.ROW_LIMIT = int(total_gaia_objects)
    # To return an unlimited number of rows set Gaia.ROW_LIMIT to -1.
    perf_start_total_query = perf_counter()
    # try:
    #     collect_antigalactic_cone(outer_cone_radius_deg=45.0)
    # except Exception:
    #     print(f"collect_antigalactic_cone failed ")
    #     pass

    try:
        collect_galactic_cone(subcone_radius_deg=90)
    except Exception:
        print(f"collect_galactic_cone failed ")
        pass

    try:
        collect_antigalactic_cone(subcone_radius_deg=90)
    except Exception:
        print(f"collect_antigalactic_cone failed ")
        pass

    # try:
    #     # anti-galactic line
    #     axis_lon, axis_lat = (180, 0)
    #     outfile_prefix="./data/antigalactic"
    #     # collect_galactic_cone(outer_cone_radius_deg=45.0)
    #     collect_outer_cone(axis_lon=axis_lon, axis_lat=axis_lat,
    #                        inner_cone_radius_deg=45.0, outer_cone_radius_deg=60.0,
    #                        outfile_prefix=outfile_prefix
    #                        )
    #     collect_outer_cone(axis_lon=axis_lon, axis_lat=axis_lat,
    #                        inner_cone_radius_deg=60.0, outer_cone_radius_deg=90.0,
    #                        outfile_prefix=outfile_prefix
    #                        )
    # except Exception:
    #     print(f"collect antigalactic failed ")
    #     pass
    #
    # try:
    #     axis_lon, axis_lat = (0, 0)
    #     outfile_prefix="./data/galactic"
    #     collect_outer_cone(axis_lon=axis_lon, axis_lat=axis_lat,
    #                        inner_cone_radius_deg=45.0, outer_cone_radius_deg=60.0,
    #                        outfile_prefix=outfile_prefix
    #                        )
    #     collect_outer_cone(axis_lon=axis_lon, axis_lat=axis_lat,
    #                        inner_cone_radius_deg=60.0, outer_cone_radius_deg=90.0,
    #                        outfile_prefix=outfile_prefix
    #                        )
    # except Exception:
    #     print(f"collect galactic failed ")
    #     pass


    print(f"total {float(Gaia.ROW_LIMIT):0.1e} query >>> elapsed: {(perf_counter() - perf_start_total_query):0.2f} seconds")


if __name__ == "__main__":
    main()
