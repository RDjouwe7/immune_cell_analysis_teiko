import streamlit as st
import pandas as pd
import sqlite3
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

elif page == "Subset Analysis":
    st.title("Baseline Melanoma Sample Analysis")
    st.markdown("""
    Melanoma PBMC samples at baseline (time=0) from patients treated with miraclib
    """)