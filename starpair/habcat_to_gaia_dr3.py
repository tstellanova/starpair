import csv
from time import perf_counter

import astropy
from astroquery.gaia import Gaia

input_file = './habcat/habcat_table.csv'  # Your input CSV file
output_file = './habcat/gaia_dr3_habcat.csv'  # Output file to save results

# Store all Hipparcos IDs and lines for later processing
hipparcos_ids = []
input_lines = []

with open(input_file, mode='r') as infile:
    csv_reader = csv.reader(infile)

    # Skip the header (if any) and store input data
    habcat_header = next(csv_reader)

    for row in csv_reader:
        hipparcos_ids.append(int(row[0]))  # Assuming 'HIP' is the first column
        input_lines.append(row)

# print(f"IDS:\n{hipparcos_ids}")

chunk_size = 1000
chunk_offset = 0
Gaia.ROW_LIMIT = 2*chunk_size

total_habcat_nodes = len(hipparcos_ids)
print(f"num HabCat nodes: {total_habcat_nodes}")
perf_start_adql_query = perf_counter()

merged_table = None

while chunk_offset < total_habcat_nodes:
    chunk_limit = min(chunk_offset + chunk_size , total_habcat_nodes)
    subset_hip_ids = hipparcos_ids[chunk_offset:chunk_limit]
    hip_ids_str = ', '.join(map(str, subset_hip_ids))

    hip_id_field = 'original_ext_source_id'
    gaia_id_field = 'SOURCE_ID'
    # ADQL query to select rows where HIP matches one of the given HIP IDs
    adql_query = f"""
        SELECT {gaia_id_field} , {hip_id_field}
        FROM gaiadr3.hipparcos2_best_neighbour
        WHERE original_ext_source_id IN ({hip_ids_str}) 
        ORDER BY {hip_id_field}
    """

    # Execute the query using Astroquery's Gaia module
    #n_query
    job = Gaia.launch_job(adql_query)
    result_table = job.get_results()
    # print(f"query results len: {len(result_table)} / {chunk_size}")
    if merged_table is None:
        merged_table = result_table
    else:
        merged_table = astropy.table.vstack([merged_table, result_table])

    # merged_table: 14192 / 17129
    print(f"merged_table: {len(merged_table)} ")
    chunk_offset = chunk_limit

print(f"merged_table: {len(merged_table)} / {total_habcat_nodes} ")

merged_table.write(output_file, format='ascii.csv', overwrite=True)
print(f"Output written to {output_file}")


print(f"ADQL query elapsed: {perf_counter() - perf_start_adql_query } secs")


