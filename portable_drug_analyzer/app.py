#!/usr/bin/env python3

import sys
import os
import logging
import pandas as pd
from Bio import Entrez
from ivablib.pubmed4125 import search_papers
from ivablib.case_report_analyzer import analyze_papers

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email for NCBI
Entrez.email = "ghhercock@gmail.com"

def main():
    """Main function."""
    
    if len(sys.argv) != 2:
        print("Usage: python app.py <drug_name>")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    
    # Create output directory if it doesn't exist
    output_dir = f"papers/{drug_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for papers
    papers = search_papers(drug_name)
    
    # Analyze papers
    analyzed_papers = analyze_papers(papers, drug_name)
    
    # Save analysis to Excel if we have results
    if not analyzed_papers.empty:
        output_file = os.path.join(output_dir, f"{drug_name}_case_reports.xlsx")
        analyzed_papers.to_excel(output_file, index=False)
        print(f"\nAnalysis saved to {output_file}")
    else:
        print("No papers could be analyzed")

if __name__ == "__main__":
    main()
