# immune_cell_analysis_teiko

A data pipeline and interactive dashboard for analyzing how the drug miraclib
affects immune cell populations in clinical trial patients.

## How To Run

### 1. Install dependencies
```bash
make setup
```

### 2. Run the pipeline
This loads the data and generates all outputs:
```bash
make pipeline
```

### 3. Launch the dashboard
```bash
make dashboard
```
The dashboard will open at `http://localhost:8501`

---

## Database Schema

The data is stored in three tables in a SQLite database.

**subjects** — one row per patient

| Column | Description |
|--------|-------------|
| subject (PK) | Unique patient ID |
| condition | Disease (e.g. melanoma) |
| age | Patient age |
| sex | Patient sex |

**samples** — one row per collected sample

| Column | Description |
|--------|-------------|
| sample (PK) | Unique sample ID |
| subject (FK) | Links to subjects table |
| project | Project ID |
| sample_type | Type of sample (e.g. PBMC) |
| treatment | Drug received |
| response | Responded to treatment (yes/no) |
| time_from_treatment_start | Days since treatment started |

**cell_counts** — one row per cell population per sample

| Column | Description |
|--------|-------------|
| id (PK) | Auto-incremented ID |
| sample (FK) | Links to samples table |
| population | Cell type (e.g. b_cell) |
| count | Raw cell count |

### Why three tables?

i picked three tables because a patient can have multiple samples taken at different timepoints. Keeping
patient info in a separate subjects table helps me prevent repeating the age, sex,
and condition for every sample

Cell counts are stored one row per population instead of five separate columns.
This just makes filtering and grouping much simpler in SQL for me.

### How does it scale?

One good thing is that no new schema or columns are needed. If there is any new projects, samples, we just need to add rows to the table.
The same goes for the cell populations.

I think indexes could be added on sample, subject, and population for faster queries

---

## Code Structure

**load_data.py**, a file used to create the database and loads the CSV

**analysis.py**, a file with three functions that perform the following:
calculates relative frequency of each cell population per sample
compares responders vs non-responders using the Mann-Whitney U test
filters baseline melanoma samples and reports demographics and averages

**dashboard.py**, a file that runs the three page Streamlit dashboard showing results from
each part of the analysis.

Also the pipeline can be run on its own to regenerate outputs without launching the dashboard.

---

## Statistical Methods

I used the Mann-Whitney U test to compare cell population frequencies between
responders and non-responders. This test was chosen because I did not want to assume normal distribution. Otherwise I would have used the T-test

---

## Dashboard

https://github.com/RDjouwe7/immune_cell_analysis_teiko/blob/main/dashboard.py