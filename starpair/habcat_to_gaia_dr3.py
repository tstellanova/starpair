import csv
from astroquery.vizier import Vizier

# Read the input CSV file without using Pandas
input_file = './habcat/habcat_table.csv'  # Your input CSV file
output_file = './habcat/gaia_dr3_habcat.csv'  # Output file to save results

# Store all Hipparcos IDs and lines for later processing
hipparcos_ids = []
input_lines = []

with open(input_file, mode='r') as infile:
    csv_reader = csv.reader(infile)

    # Skip the header (if any) and store input data
    header = next(csv_reader)

    for row in csv_reader:
        hipparcos_ids.append(row[0])  # Assuming 'HIP' is the first column
        input_lines.append(row)

# Query Vizier for Gaia DR3 identifiers based on Hipparcos catalog numbers
vizier = Vizier(columns=['HIP', 'Source'], row_limit=-1)
catalog = "I/350/gaiaedr3"  # Gaia DR3 catalog

# Query the catalog for the Hipparcos to Gaia DR3 mapping
result = vizier.query_constraints(catalog=catalog, HIP=hipparcos_ids)

# Prepare a dictionary for fast look-up of Gaia DR3 source by HIP number
hip_to_gaia = {}
if len(result) > 0:
    for row in result[0]:
        hip_to_gaia[str(row['HIP'])] = str(row['Source'])

# Write the output CSV with Gaia DR3 source IDs and the original data
with open(output_file, mode='w', newline='') as outfile:
    csv_writer = csv.writer(outfile)

    # Write the header, adding Gaia Source ID as the first column
    csv_writer.writerow(['Gaia_DR3_Source'] + header)

    for line in input_lines:
        hip_id = line[0]
        gaia_id = hip_to_gaia.get(hip_id, 'N/A')  # Use 'N/A' if no match is found
        csv_writer.writerow([gaia_id] + line)

print(f"Output written to {output_file}")