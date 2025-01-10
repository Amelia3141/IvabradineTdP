# Ivabradine TdP Risk Assessment Dashboard

A Streamlit-based dashboard for analyzing Torsades de Pointes (TdP) risk associated with Ivabradine and other drugs.

## Features

- Risk category analysis using CredibleMeds data
- Literature review from PubMed (case reports, cohort studies, clinical trials)
- hERG channel activity analysis
- Interactive visualizations

## Installation

```bash
pip install -r requirements.txt
```

## Running Locally

```bash
streamlit run app.py
```

## Deployment

This app is deployed on Streamlit Cloud. Visit [https://ivabradine-tdp.streamlit.app](https://ivabradine-tdp.streamlit.app) to use the live version.

### Deploy Your Own Instance

1. Fork this repository
2. Go to [Streamlit Cloud](https://streamlit.io/cloud)
3. Click "New app"
4. Select your forked repository
5. Select `app.py` as the main file
6. Click "Deploy"

## Data Sources

- CredibleMeds database for TdP risk categorization
- PubMed for literature analysis
- PubChem for hERG activity data

## License

MIT License
