import pandas as pd
import sqlite3
import os


CSV_PATH = 'data/cell-count.csv'
DB_PATH = 'immune_cells.db'

def init_db(conn):
    """
       function to initialize my SQLite database
    Args:
        conn (sqlite3.Connection): Active SQLite db connection
    Returns:
        None
    """    
    controller = conn.cursor()

    #creating my 3 tables based on the features in the data

    #subjects table
    controller.execute(
        """
        CREATE TABLE IF NOT EXISTS subjects (
            subject     TEXT PRIMARY KEY,
            condition   TEXT,
            age         INTEGER,
            sex         TEXT
        )
    """
    )

    # samples table
    controller.execute(
        """
        CREATE TABLE IF NOT EXISTS samples (
            sample                 TEXT PRIMARY KEY,
            subject                     TEXT,
            project                     TEXT,
            sample_type                 TEXT,
            treatment                   TEXT,
            response                    TEXT,
            time_from_treatment_start   INTEGER,
            FOREIGN KEY (subject) REFERENCES subjects(subject)
        )
    """
    )

    # cell_counts table
    controller.execute(
        """
        CREATE TABLE IF NOT EXISTS cell_counts (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            sample   TEXT,
            population  TEXT,
            count       INTEGER,
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
        )
    """)

    conn.commit()
    print("Database tables created.")


def load_data(conn):
    """ function to load data into my database

    Args:
        conn (sqlite3.Connection): Active SQLite db connection

    Returns:
        None
    """    
    df = pd.read_csv(CSV_PATH)

     # loading subjects and dropping duplicates
    subjects_df = df[["subject", "condition", "age", "sex"]].drop_duplicates()
    subjects_df.to_sql("subjects", conn, if_exists="replace", index=False)
    print("Loaded subjects")

    # loading the samples
    samples_df = df[[
        "sample", "subject", "project", "sample_type",
        "treatment", "response", "time_from_treatment_start"
    ]].copy()
    samples_df.to_sql("samples", conn, if_exists="replace", index=False)
    print(f"Loaded {len(samples_df)} samples")


     # loading cell_counts
    cell_cols = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    counts_df = df[["sample"] + cell_cols].copy()

    # melting the df to make my analysis easier
    counts_long = counts_df.melt(
        id_vars="sample",
        value_vars=cell_cols,
        var_name="population",
        value_name="count"
    )
    counts_long.to_sql("cell_counts", conn, if_exists="replace", index=False)
    print("Loaded cell_counts records")



if __name__ == "__main__":

    conn = sqlite3.connect(DB_PATH)
    init_db(conn)
    load_data(conn)
    conn.close()
    print("Data saved to", DB_PATH)