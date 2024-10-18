import csv
import os
import sqlite3
import numpy as np
import argparse


g_star_by_gaia_id = {}

def infer_type(value):
    """Helper function to infer the SQLite type of a CSV value."""
    try:
        int(value)
        return "INTEGER"
    except ValueError:
        try:
            float(value)
            return "REAL"
        except ValueError:
            return "TEXT"

def create_table_from_csv(db_name, csv_file, table_name, primary_key):
    """Create a table in an SQLite database from a CSV file and import the data."""
    # Connect to SQLite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    input_header = None
    field_types = None
    # Scan the first few lines of the CSV file to infer field types
    with open(csv_file, 'r') as input_file:
        reader = csv.reader(input_file)
        input_header = next(reader)

        # Read the first row to infer the types
        first_row = next(reader)
        field_types = [infer_type(value) for value in first_row]
        print(f"field_types: {field_types}")

        # Construct the SQL statement to create the table
        fields = [
            f'"{col}" {ftype}' + (f' PRIMARY KEY' if col == primary_key else '')
            for col, ftype in zip(input_header, field_types)
        ]
        create_table_sql = f'CREATE TABLE IF NOT EXISTS "{table_name}" ({", ".join(fields)});'

        # Create the table
        print(f"Create table sql: {create_table_sql}")
        cursor.execute(create_table_sql)

# CREATE TABLE IF NOT EXISTS "thzc" ("TIC_ID" INTEGER, "Teff" INTEGER, "Tmag" REAL, "r_est" REAL, "Obs" REAL, "a_RV" REAL, "a_EA" REAL, "a_EM" REAL, "Per_RV" REAL, "Per_EA" REAL, "Per_EM" REAL, "Elat" REAL, "Glat" REAL, "Gaia_ID" INTEGER PRIMARY KEY, "RA" REAL, "Dec" REAL, "asep_RV" REAL, "asep_EA" REAL, "asep_EM" REAL, "Ttime" REAL, "TWOMASS" TEXT, "HZ_flag" INTEGER, "JWST_flag" INTEGER, "Earth_flag" INTEGER, "S1_flag" INTEGER, "S2_flag" INTEGER, "S3_flag" INTEGER, "S4_flag" INTEGER, "S5_flag" INTEGER, "S6_flag" INTEGER, "" TEXT);


# Reopen the CSV file to import the data
    with open(csv_file, 'r') as input_file:
        reader = csv.DictReader(input_file)

        # Construct the SQL statement to insert the data
        placeholders = ', '.join('?' for _ in input_header)
        insert_query = f'INSERT INTO "{table_name}" ({", ".join(input_header)}) VALUES ({placeholders})'

        # Insert all rows into the database
        for row in reader:
            # TODO use field_types for speed?
            gaia_id = np.uint64(row[primary_key])
            existing_gaia = g_star_by_gaia_id.get(gaia_id)
            if existing_gaia is None:
                g_star_by_gaia_id[gaia_id] = gaia_id
                values = [row[col] if infer_type(row[col]) == 'TEXT' else (int(row[col]) if infer_type(row[col]) == 'INTEGER' else float(row[col])) for col in input_header]
                cursor.execute(insert_query, values)
            else:
                print(f"ignoring existing Gaia ID: {gaia_id}")


    # Commit
    conn.commit()
    cursor = None
    return conn



# def open_db(db_file_path: str):
#     '''
#     Connect to SQLite (creates a new database if it doesn't exist)
#     :param db_file_path: Path to the db file
#     :return:  db connection, cursor
#     '''
#     conn = sqlite3.connect(db_file_path)
#     cursor = conn.cursor()
#
#     # Create a table with the specified fields and composite primary key (src_id, dst_id)
#     cursor.execute('''
#         CREATE TABLE IF NOT EXISTS neighbor_pairs (
#             max_node_dist REAL,
#             min_node_dist REAL,
#             ang_sep REAL,
#             radial_sep REAL,
#             src_coord_lon REAL,
#             src_coord_lat REAL,
#             dst_coord_lon REAL,
#             dst_coord_lat REAL,
#             src_id_str TEXT,
#             dst_id_str TEXT,
#             src_id INTEGER,
#             dst_id INTEGER,
#             PRIMARY KEY (src_id, dst_id)
#         )
#     ''')
#
#     return conn, cursor

def main():
    parser = argparse.ArgumentParser(description='Import THZC catalog into a sqlite db')
    parser.add_argument('src_path', nargs='?',
                        default="./tess/tess_hab_zone_cat_all.fits.csv",
                        # default="./tess/tess_hab_zone_cat_d100.fits.csv",
                        help="TESS Hab Zone Catalog file with `.csv` extension",
                        )
    parser.add_argument('-o', dest='outpath', type=str,
                        help='Output path for db')

    args = parser.parse_args()
    data_file_name = args.src_path
    basename_sans_ext = os.path.splitext(os.path.basename(data_file_name))[0]

    db_file_path = args.outpath
    if db_file_path is None:
        db_file_path = f"./tess/{basename_sans_ext}_raw.db"
    print(f"Loading: {data_file_name} , output to: {db_file_path}")
    conn = create_table_from_csv(db_file_path, data_file_name, "thzc", 'Gaia_ID')
    cursor = conn.cursor()
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_r_est ON thzc (r_est)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_ra ON thzc (RA)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_dec ON thzc (Dec)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_tic_id ON thzc (TIC_ID)')
    conn.commit()
    conn.close()


    print(f"CSV data imported into the SQLite database at: \n{db_file_path}")

if __name__ == "__main__":
    main()
