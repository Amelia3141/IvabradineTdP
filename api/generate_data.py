import sys
import os
import json
from datetime import datetime
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ivablib.case_report_analyzer import CaseReportAnalyzer
from ivablib.pubmed4125 import PubMedSearcher

def generate_drug_data(drug_name):
    """Generate analysis data for a given drug"""
    searcher = PubMedSearcher()
    analyzer = CaseReportAnalyzer()
    
    # Search for papers
    papers = searcher.search_drug(drug_name)
    
    # Analyze papers
    results = analyzer.analyze_papers(papers)
    
    # Convert to JSON-serializable format
    output = {
        "drug_name": drug_name,
        "timestamp": datetime.utcnow().isoformat(),
        "paper_count": len(papers),
        "analysis": results,
    }
    
    # Save to api/data directory
    os.makedirs("api/data", exist_ok=True)
    output_file = f"api/data/{drug_name.lower()}.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)
        
    # Update drug list
    drug_list_file = "api/data/drugs.json"
    try:
        with open(drug_list_file, "r") as f:
            drugs = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        drugs = {"drugs": []}
        
    if drug_name.lower() not in [d.lower() for d in drugs["drugs"]]:
        drugs["drugs"].append(drug_name)
        with open(drug_list_file, "w") as f:
            json.dump(drugs, f, indent=2)
            
    return output

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_data.py <drug_name>")
        sys.exit(1)
    generate_drug_data(sys.argv[1])
