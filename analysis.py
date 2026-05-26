import sqlite3
import pandas as pd
import os

DB_PATH = "immune_cells.db"
OUTPUT_DIR = "outputs"


def get_connection():
    """ function to connect to the SQLite database

    Returns:
        sqlite3.Connection: Return a connection to the SQLite database
    """    
    return sqlite3.connect(DB_PATH)


def frequency_table(conn):
    """ Calculate relative frequency of each cell population per sample.
    
    Args:
        conn (sqlite3.Connection): Active SQLite db connection

    Returns:
        freq_table (dataframe): Returns a dataframe with sample, total_count, population, count, percentage
    """    
    query = """
        SELECT
            A.sample,
            A.subject,
            A.treatment,
            A.response,
            A.sample_type,
            A.time_from_treatment_start,
            C.condition,
            C.sex,
            B.population,
            B.count
        FROM samples A
        JOIN cell_counts B ON A.sample = B.sample
        JOIN subjects C ON A.subject = C.subject
    """
    df = pd.read_sql_query(query, conn)

    # calculating total count per sample
    total = df.groupby("sample")["count"].sum().reset_index()
    total = total.rename(columns={"count": "total_count"})

    # merging total back into main df
    df = df.merge(total, on="sample")

    # calculating the percentages
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(2)

    freq_table = df[[
        "sample", "total_count", "population", "count", "percentage"
    ]]

    return freq_table


if __name__ == "__main__":

    #creating the output folder if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    conn = get_connection()
    
    summary_table = frequency_table(conn)

    #putting my freq table in csv
    summary_table.to_csv(f"{OUTPUT_DIR}/summary_table.csv", index=False)
    print("Summary table saved")
    conn.close()