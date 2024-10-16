"""
Import starpairs into a database,
culling pairs where neither node appears
in a list of habitable star systems

"""
import csv
import os
import sqlite3
import numpy as np
import argparse

def open_db(db_file_path: str):
    '''
    Connect to SQLite (creates a new database if it doesn't exist)
    :param db_file_path: Path to the db file
    :return:  db connection, cursor
    '''
    conn = sqlite3.connect(db_file_path)
    cursor = conn.cursor()

    # Create a table with the specified fields and composite primary key (src_id, dst_id)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS neighbor_pairs (
            max_node_dist REAL,
            min_node_dist REAL,
            ang_sep REAL,
            radial_sep REAL,
            src_coord_lon REAL,
            src_coord_lat REAL,
            dst_coord_lon REAL,
            dst_coord_lat REAL,
            src_id_str TEXT,
            dst_id_str TEXT,
            src_id INTEGER,
            dst_id INTEGER,
            PRIMARY KEY (src_id, dst_id)
        )
    ''')

    return conn, cursor

# Function to split coordinate string and convert to float64 values
def parse_coordinate_str(coord_str)->(np.float64, np.float64):
    lon, lat = coord_str.split()
    return np.float64(float(lon)), np.float64(float(lat))

def main():
    parser = argparse.ArgumentParser(description='Convert star pairs to db')
    parser.add_argument('src_path', nargs='?',
                        # default="./data/galactic_L2e9_r90_d100_0_0_ma250_spairs.csv",
                        default="./data/antigalactic_L2e9_r90_d100_18000_0_ma250_spairs.csv",
                        help="Source data file with `.csv` extension",
                        )
    parser.add_argument('-f', dest='hab_path',
                        default="./tess/tic_to_gaia_mapping.csv",
                        help="csv habitability catalog file with rows: external ID, gaia source ID",
                        )
    parser.add_argument('-o', dest='outpath', type=str,
                        help='Output path for db')

    args = parser.parse_args()
    data_file_name = args.src_path
    basename_sans_ext = os.path.splitext(os.path.basename(data_file_name))[0]
    filter_file_path = args.hab_path

    # Initialize an empty dictionary
    habitable_map = {}
    hab_list_row_count = 0
    # import the habitable star system list file
    with open(filter_file_path, mode='r') as hab_list_file:
        csv_reader = csv.reader(hab_list_file)
        unused_header = next(csv_reader)
        for hab_list_row in csv_reader:
            hab_list_row_count += 1
            hab_list_uid, gaia_dr3_id, = hab_list_row
            habitable_map[gaia_dr3_id] = hab_list_uid
        hab_list_file.close()
        hab_list_file = None
        hab_list_row = None
    print(f"num habitable entries: {hab_list_row_count} vs {len(habitable_map)}")

    db_file_path = args.outpath
    if db_file_path is None:
        db_file_path = f"./data/{basename_sans_ext}_hab_filt.db"
    print(f"Input: {data_file_name} \nOutput: {db_file_path}")

    conn, cursor = open_db(db_file_path)

    uninhab_skip_count = 0
    kept_pairs_count = 0
    # Read CSV and insert data into the table
    with open(data_file_name, 'r') as in_file:
        csv_reader = csv.reader(in_file)
        headers = next(csv_reader)  # Skip the first header row
        for row in csv_reader:
            max_node_dist = np.float64(row[0])
            min_node_dist = np.float64(row[1])
            ang_sep = np.float64(row[2])
            radial_sep = np.float64(row[3])

            # Split source and destination coordinates
            src_coord_lon, src_coord_lat = parse_coordinate_str(row[4])
            dst_coord_lon, dst_coord_lat = parse_coordinate_str(row[5])

            src_id_str = row[6]
            dst_id_str = row[7]

            # cull star pairs that are not in the habitable list
            if habitable_map.get(src_id_str) is None and habitable_map.get(dst_id_str) is None:
                uninhab_skip_count += 1
                continue
            src_id = np.uint64(src_id_str)
            dst_id = np.uint64(dst_id_str)

            kept_pairs_count += 1
            # Insert the data into the table
            cursor.execute('''
                INSERT OR IGNORE INTO neighbor_pairs (
                    max_node_dist, min_node_dist, ang_sep, radial_sep,
                    src_coord_lon, src_coord_lat, dst_coord_lon, dst_coord_lat,
                    src_id_str, dst_id_str, src_id, dst_id
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (max_node_dist, min_node_dist, ang_sep, radial_sep,
                  src_coord_lon, src_coord_lat, dst_coord_lon, dst_coord_lat,
                  src_id_str, dst_id_str, src_id, dst_id))

    # Commit the changes to the database
    conn.commit()

    # Create indexes for faster lookups on src_id and dst_id
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_src_id ON neighbor_pairs (src_id)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_dst_id ON neighbor_pairs (dst_id)')

    # Close the connection
    conn.close()

    print(f"total input star pairs: {kept_pairs_count + uninhab_skip_count}")
    print(f"uninhab_skip_count {uninhab_skip_count} ")
    print(f"kept_pairs_count: {kept_pairs_count}")
    print(f"CSV data imported into the SQLite database at: \n{db_file_path}")

if __name__ == "__main__":
    main()
