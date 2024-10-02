import glob
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np

filepath_stem = "./data/antigalactic_L250_r10_d306"
# Define the filename pattern (e.g., all .fits files in the current directory)
fits_files = glob.glob(f"{filepath_stem}*.fits")  # You can modify the pattern as needed
print(f"fits_files: {fits_files}")

# Initialize an empty list to store tables
tables = []

# Loop over the FITS files to read the tables
for filename in fits_files:
    with fits.open(filename) as hdul:
        # Assuming the table is in the first HDU
        table = Table(hdul[1].data)
        tables.append(table)

# Merge tables using vstack
merged_table = vstack(tables)

# Remove duplicates based on the 'DESIGNATION' column
unique_table = merged_table.group_by('DESIGNATION')
unique_table = unique_table.groups.aggregate(np.unique)

# # Output the result
# print(unique_table)

# Save the merged unique table back to a new FITS file if needed
unique_table.write(f"{filepath_stem}_merged.fits", overwrite=True)