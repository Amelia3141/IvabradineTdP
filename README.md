# Ivabradine TdP Risk Analyzer

A Streamlit web application for analyzing the Torsades de Pointes (TdP) risk of Ivabradine through hERG channel analysis and literature review.

## Research Context

This project emerged from my undergraduate dissertation research on Ivabradine and its potential arrhythmogenic risks. The presentation below provides the academic foundation and findings that motivated the development of this analytical tool.

### 📊 Thesis Presentation: "Ivabradine and the Risks of Torsades: a Case-Based Analysis"

**[📄 View Full Presentation (PDF)](Ivabradine%20Thesis%20Presentation.pdf)**

*Supervisor: Prof. J Hancox | University of Bristol | 2024*

**Key Research Findings:**
- Systematic analysis of Ivabradine's cardiac safety profile
- Investigation of TdP risk factors and mechanisms
- Development of risk assessment methodologies
- Literature review of case reports and clinical data

This web application translates the research findings into a practical tool for automated TdP risk assessment, incorporating the methodologies developed during the thesis work.

---

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

## Academic Background

This project builds upon comprehensive thesis research that:
- Analyzed clinical case reports of Ivabradine-associated arrhythmias
- Developed BioBERT-based text mining systems for pharmacovigilance
- Created automated systematic review tools for torsadogenic medications
- Achieved a dissertation grade of 84 overall with 89 on Viva examination

The web application implements the risk assessment frameworks developed during this research, making them accessible for clinical and research applications.

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

---

## Citation

If you use this tool in your research, please cite:

```
Ghanea, A. (2024). Ivabradine and the Risk of Torsades de Pointes: A Case-Based Analysis. 
University of Bristol. [Thesis]
```
