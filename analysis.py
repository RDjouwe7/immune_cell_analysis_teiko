import sqlite3
import pandas as pd
import os
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

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
        "sample", "total_count", "population", "count", "percentage",
        "subject", "treatment", "response", "sample_type",
        "time_from_treatment_start", "condition", "sex"
    ]]

    return freq_table

def statistical_analysis(table):
    """
    Compare cell population frequencies between responders and non-responders.
    This analysis only includes melanoma patients treated with miraclib using PBMC samples

    Args:
        table (pandas.DataFrame): Dataframe containing frequency analysis results and sample metadata
    Returns:
        results_df (pandas.DataFrame): Statistical test results for each immune cell population
    """    
    # filtering to get the desired data
    filtered= table[
        (table["condition"] == "melanoma") &
        (table["treatment"] == "miraclib") &
        (table["sample_type"] == "PBMC") &
        (table["response"].isin(["yes", "no"]))
    ].copy()

    # creating my boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        data=filtered,
        x="population",
        y="percentage",
        hue="response",
        palette={"yes": "green", "no": "red"}
    )
    plt.title("Cell Population Frequencies: Responders vs Non-Responders\n"
              "(Melanoma patients on miraclib, PBMC samples)")
    plt.xlabel("Cell Population")
    plt.ylabel("Relative Frequency (%)")
    plt.legend(title="Response")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/results.png", dpi=150)
    plt.close()
    print("Boxplot done")

    # doing my statistical testing
    results = []
    for population in filtered["population"].unique():
        pop_data = filtered[filtered["population"] == population]

        responders = pop_data[pop_data["response"] == "yes"]["percentage"]
        non_responders = pop_data[pop_data["response"] == "no"]["percentage"]

        # performing mann-whitney u test in case of not normal distribution
        stat, pvalue = stats.mannwhitneyu(
            responders,
            non_responders,
            alternative="two-sided" 
        )

        results.append({
            "population": population,
            "mean_pct_responders": round(responders.mean(), 2),
            "mean_pct_non_responders": round(non_responders.mean(), 2),
            "p_value": round(pvalue, 4),
            "significant": "yes" if pvalue < 0.05 else "no"
        })

    results_df = pd.DataFrame(results).sort_values("p_value")
    results_df.to_csv(f"{OUTPUT_DIR}/statistical_results.csv", index=False)

    return results_df

def subset_analysis(conn):
    """"
    function to analyze baseline melanoma PBMC samples treated with miraclib.

    Args:
        conn (sqlite3.Connection): Active SQLite database connection.
        
    Returns:
        tuple:
            - sample_per_project (pandas.DataFrame): Sample counts per project
            - response_counts (pandas.DataFrame): Subject counts by response
            - sex_counts (pandas.DataFrame): Subject counts by sex
    """    
    # base query needs melanoma, PBMC, miraclib, time=0
    base_query = """
        SELECT 
            A.sample,
            A.subject,
            A.project,
            A.treatment,
            A.response,
            A.sample_type,
            A.time_from_treatment_start,
            B.condition,
            B.sex
        FROM samples A
        JOIN subjects B ON A.subject = B.subject
        WHERE B.condition = 'melanoma'
        AND A.sample_type = 'PBMC'
        AND A.treatment = 'miraclib'
        AND A.time_from_treatment_start = 0
    """
    df = pd.read_sql_query(base_query, conn)

    #Total number of baseline samples (melanoma, PBMC, miraclib)
    print(f"\nTotal baseline samples found: {len(df)}")

    # getting my samples per project 
    samples_per_project = df.groupby("project")["sample"].count().reset_index()
    print(samples_per_project)
    samples_per_project.columns = ["project", "sample_count"]
    print("\nSamples per project:")
    print(samples_per_project.to_string(index=False))

    # getting subjects per response group
    # using nunique to avoid counting same subject twice
    response_counts = df.groupby("response")["subject"].nunique().reset_index()
    response_counts.columns = ["response", "subject_count"]
    print("\nSubjects per response:")
    print(response_counts.to_string(index=False))

    # getting the subjects per sex
    sex_counts = df.groupby("sex")["subject"].nunique().reset_index()
    sex_counts.columns = ["sex", "subject_count"]
    print("\nSubjects per sex:")
    print(sex_counts.to_string(index=False))


    samples_per_project.to_csv(f"{OUTPUT_DIR}/samples_per_project.csv", index=False)
    response_counts.to_csv(f"{OUTPUT_DIR}/response_counts.csv", index=False)
    sex_counts.to_csv(f"{OUTPUT_DIR}/sex_counts.csv", index=False)


    return samples_per_project, response_counts, sex_counts

if __name__ == "__main__":

    #creating the output folder if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    conn = get_connection()
    
    summary_table = frequency_table(conn)

    #putting my freq table in csv
    summary_table.to_csv(f"{OUTPUT_DIR}/summary_table.csv", index=False)
    print("Summary table saved")

    results = statistical_analysis(summary_table)
    print("\nStatistical Results:")
    print(results.to_string(index=False))
    samples_per_project, response_counts, sex_counts = subset_analysis(conn)

    conn.close()