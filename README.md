# Ivabradine TdP Risk Analyzer

A Streamlit web application for analyzing the Torsades de Pointes (TdP) risk of Ivabradine through hERG channel analysis and literature review.

## Features

- **hERG Channel Analysis**: 
  - Retrieves hERG IC50 values from PubChem bioassays and bioactivity data
  - Calculates theoretical and plasma concentrations
  - Estimates TdP risk based on concentration ratios
- **Literature Search**:
  - Searches PubChem for ion channel data
  - Analyzes bioassay results
  - Retrieves bioactivity data
- **Concentration Analysis**:
  - Calculates theoretical maximum concentrations
  - Estimates plasma concentrations with 40% bioavailability
  - Computes hERG IC50 ratios

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Amelia3141/IvabradineTdP.git
cd IvabradineTdP
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Configure NCBI credentials:
Create a `.streamlit/secrets.toml` file with your NCBI credentials:
```toml
NCBI_EMAIL = "your.email@example.com"
NCBI_API_KEY = "your_api_key"
```

4. Run the application:
```bash
streamlit run app.py
```

## Usage

1. Enter the drug name (e.g., "ivabradine")
2. Input the doses to analyze
3. View the analysis results:
   - hERG IC50 value and source
   - Concentration calculations
   - TdP risk assessment
   - Citations and references

## Dependencies

- streamlit
- plotly
- pandas
- requests
- beautifulsoup4
- PyPDF2

## Recent Updates

- Switched to PubChem API for more comprehensive hERG data
- Added bioassay and bioactivity searches
- Improved concentration calculations
- Enhanced error handling and logging

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
