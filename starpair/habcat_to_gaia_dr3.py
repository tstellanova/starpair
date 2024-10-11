import csv
from time import perf_counter

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

n_query = len(hipparcos_ids)
print(f"num hipparcos_ids: {n_query}")
subset_hip_ids = hipparcos_ids
hip_ids_str = ', '.join(map(str, subset_hip_ids))

hip_id_field = 'original_ext_source_id'
gaia_id_field = 'SOURCE_ID'
# ADQL query to select rows where HIP matches one of the given HIP IDs
# Note we order by HIP, which is the same order as the original habcat table,
# allowing for easy verification
adql_query = f"""
    SELECT {hip_id_field}, {gaia_id_field} 
    FROM gaiadr3.hipparcos2_best_neighbour
    WHERE original_ext_source_id IN ({hip_ids_str}) 
    ORDER BY {hip_id_field}
"""
#     ORDER BY HIP

# Execute the query using Astroquery's Gaia module
perf_start_adql_query = perf_counter()
Gaia.ROW_LIMIT = n_query
job = Gaia.launch_job(adql_query)
result_table = job.get_results()
result_table.write(output_file, format='ascii.csv', overwrite=True)

# Extract the HIP and corresponding Gaia source IDs from the results
result_rows = zip(result_table[hip_id_field], result_table[gaia_id_field])
print(f"ADQL query elapsed: {perf_counter() - perf_start_adql_query } secs")
hip_to_gaia_dict = dict(result_rows)


print(f"Output written to {output_file}")