# TdP Risk Predictor

A tool for analyzing drug-induced Torsades de Pointes (TdP) risk by extracting and analyzing case reports from medical literature.

## Features

- Automated PubMed search for drug case reports
- PDF download and text extraction
- Analysis of key clinical parameters:
  - Patient age and sex
  - Drug dosage
  - Heart rate
  - Blood pressure
  - QTc measurements
- Excel report generation with findings

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Amelia3141/TdPRiskPredictor.git
cd TdPRiskPredictor
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the analyzer with a drug name:
```bash
python app.py donepezil
```

The tool will:
1. Search PubMed for case reports
2. Download and analyze the papers
3. Extract key clinical information
4. Generate an Excel report with the findings

## Project Structure

```
TdPRiskPredictor/
├── app.py              # Main entry point
├── requirements.txt    # Python dependencies
├── data/
│   └── drugnames.csv  # Drug name mappings
└── ivablib/           # Core library
    ├── __init__.py
    ├── case_report_analyzer.py  # Text analysis
    └── pubmed4125.py           # PubMed interface
```

## Dependencies

- pandas
- numpy
- requests
- beautifulsoup4
- PyPDF2
- biopython

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
