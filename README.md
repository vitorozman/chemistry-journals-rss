# chemistry-journals-rss

A staging repository for exploring and visualising chemistry research activity in Texas. It fetches the latest publications from top chemistry journals via RSS/PubMed feeds, identifies potential Texas-affiliated researcher candidates, and surfaces statistics useful for further investigation.

> **Status:** staging / prototype — used for idea validation and exploratory analysis.

---

## Live Apps

| App | Description | Link |
|---|---|---|
| **What's New** | Latest publications from top chemistry journals (last 2 weeks), grouped by journal | [whats-new-texas-chemistry.streamlit.app](https://whats-new-texas-chemistry.streamlit.app) |
| **Texas Researchers** | Candidate researchers with Texas affiliation not yet in the database, with their publications | [texas-researchers.streamlit.app](https://texas-researchers.streamlit.app) |

---


**Texas Researchers** candidates satisfy all of:
- Have a **Texas affiliation**
- Are **not** currently in the database
- Have ≥ 1 publication in a **top-rated chemistry journal** where they appear as **first or last author** with their **full name** listed

---

## Stack

- **Python** — data fetching, parsing, enrichment
- **Streamlit** — interactive web UI (no frontend code required)
- **Streamlit Cloud** — hosting for both live apps
