import streamlit as st
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
import os

# page config
st.set_page_config(
    page_title="Immune Cell Analysis",
    layout="wide"
)

DB_PATH = "immune_cells.db"
OUTPUT_DIR = "outputs"


def get_connection():
    """ function to connect to the SQLite database

    Returns:
        sqlite3.Connection: Return a connection to the SQLite database
    """    
    return sqlite3.connect(DB_PATH)


# setting my sidebar navigation
st.sidebar.title("Immune Cell Analysis")
st.sidebar.markdown("Loblaw Bio Clinical Trial Dashboard")
page = st.sidebar.radio(
    "Navigate",
    ["Overview", "Statistical Analysis", "Subset Analysis"]
)

if page == "Overview":
    st.title("Cell Population Frequencies")
    st.markdown("Relative frequency of each immune cell population per sample")

    summary = pd.read_csv(f"{OUTPUT_DIR}/summary_table.csv")

    # my metrics
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Samples", summary["sample"].nunique())
    col2.metric("Total Subjects", summary["subject"].nunique())
    col3.metric("Cell Populations", summary["population"].nunique())

    # the table
    st.dataframe(
        summary[["sample", "population", "count", "percentage", "total_count"]],
        use_container_width=True
    )

    # adding a download button
    st.download_button(
        label="Download Summary Table",
        data=summary.to_csv(index=False),
        file_name="summary_table.csv",
        mime="text/csv"
    )

elif page == "Statistical Analysis":
    st.title("Responders vs Non-Responders")
    st.markdown("""
    Comparing cell population frequencies between responders and non-responders
    in melanoma patients treated with miraclib (PBMC samples only)
    """)

    results = pd.read_csv(f"{OUTPUT_DIR}/statistical_results.csv")
    summary = pd.read_csv(f"{OUTPUT_DIR}/summary_table.csv")

    # showing the boxplot
    st.subheader("Cell Population Frequencies")
    filtered = summary[
        (summary["condition"] == "melanoma") &
        (summary["treatment"] == "miraclib") &
        (summary["sample_type"] == "PBMC") &
        (summary["response"].isin(["yes", "no"]))
    ]

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(
        data=filtered,
        x="population",
        y="percentage",
        hue="response",
        palette={"yes": "green", "no": "red"},
        ax=ax
    )
    ax.set_title("Cell Population Frequencies: Responders vs Non-Responders")
    ax.set_xlabel("Cell Population")
    ax.set_ylabel("Relative Frequency (%)")
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()

    # showing statistical results
    st.subheader("Statistical Results Using Mann-Whitney U Test")
    st.markdown("P-value < 0.05 indicates a statistically significant difference.")

    # highlighting significant rows
    def highlight_significant(row):
        """
        Highlight rows with statistically significant results.
        Args:
            row (pandas.Series): Row from the statistical results dataframe.
        Returns:
            list: List of CSS styles for dataframe row formatting
        """        
        if row["significant"] == "yes":
            return ["background-color: lightgreen"] * len(row)
        return [""] * len(row)

    styled = results.style.apply(highlight_significant, axis=1)
    st.dataframe(styled, use_container_width=True)

elif page == "Subset Analysis":
    st.title("Baseline Melanoma Sample Analysis")
    st.markdown("""
    Melanoma PBMC samples at baseline (time=0) from patients treated with miraclib
    """)
    conn = get_connection()

    # samples per project
    st.subheader("Samples Per Project")
    samples_per_project = pd.read_csv(f"{OUTPUT_DIR}/samples_per_project.csv")
    st.dataframe(samples_per_project, use_container_width=True)

    # response and sex counts side by side
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Subjects by Response")
        response_counts = pd.read_csv(f"{OUTPUT_DIR}/response_counts.csv")
        st.dataframe(response_counts, use_container_width=True)

    with col2:
        st.subheader("Subjects by Sex")
        sex_counts = pd.read_csv(f"{OUTPUT_DIR}/sex_counts.csv")
        st.dataframe(sex_counts, use_container_width=True)

    conn.close()