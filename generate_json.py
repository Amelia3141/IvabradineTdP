import json
import os
import pandas as pd
from ivablib.case_report_analyzer import CaseReportAnalyzer
from ivablib.pubmed4125 import search_papers

def analyze_drug(drug_name):
    # Initialize analyzer
    analyzer = CaseReportAnalyzer()
    
    # Fetch papers for the drug
    papers = search_papers(drug_name)
    
    # Analyze papers
    case_reports_df = analyzer.analyze_papers(papers, drug_name)
    
    # Convert DataFrame to desired JSON structure
    case_reports = []
    for _, row in case_reports_df.iterrows():
        report = {
            "title": row.get("title", ""),
            "authors": row.get("authors", ""),
            "year": int(row.get("year", 0)) if pd.notna(row.get("year")) else None,
            "journal": row.get("journal", ""),
            "doi": row.get("doi", ""),
            "age": float(row.get("age")) if pd.notna(row.get("age")) else None,
            "sex": row.get("sex"),
            "qtc": float(row.get("qtc")) if pd.notna(row.get("qtc")) else None,
            "qt_uncorrected": float(row.get("qt_uncorrected")) if pd.notna(row.get("qt_uncorrected")) else None,
            "heart_rate": float(row.get("heart_rate")) if pd.notna(row.get("heart_rate")) else None,
            "blood_pressure": row.get("blood_pressure"),
            "medical_history": row.get("medical_history"),
            "medication_history": row.get("medication_history"),
            "treatment_course": row.get("treatment_course"),
            "outcome": "TdP" if row.get("had_tdp") == "Yes" else "Unknown"
        }
        case_reports.append(report)
    
    # Calculate statistics
    ages = [r["age"] for r in case_reports if r["age"] is not None]
    age_stats = {
        "mean": sum(ages) / len(ages) if ages else 0,
        "std": 0,  # TODO: Calculate standard deviation
        "min": min(ages) if ages else 0,
        "max": max(ages) if ages else 0
    }
    
    sex_counts = {"male": 0, "female": 0}
    for report in case_reports:
        if report["sex"] and report["sex"].lower() in ["male", "female"]:
            sex_counts[report["sex"].lower()] += 1
    
    # Count risk factors (using medical history as risk factors)
    risk_factors = {}
    for _, row in case_reports_df.iterrows():
        if pd.notna(row.get("medical_history")):
            conditions = [c.strip() for c in str(row["medical_history"]).split(",")]
            for condition in conditions:
                if condition:
                    risk_factors[condition] = risk_factors.get(condition, 0) + 1
    
    analysis = {
        "drug_name": drug_name,
        "paper_count": len(case_reports),
        "analysis": {
            "age_stats": age_stats,
            "sex_distribution": sex_counts,
            "risk_factors": risk_factors,
            "case_reports": case_reports
        }
    }
    
    return analysis

def main():
    # Create data directory if it doesn't exist
    os.makedirs("data", exist_ok=True)
    
    # List of drugs to analyze
    drugs = ["donepezil", "ivabradine", "diphenhydramine", "azithromycin", "ranolazine", "aripiprazole"]
    
    # Generate JSON for each drug
    for drug in drugs:
        print(f"Analyzing {drug}...")
        analysis = analyze_drug(drug)
        
        with open(f"data/{drug}.json", "w") as f:
            json.dump(analysis, f, indent=2)
    
    # Create drugs.json with list of available drugs
    with open("data/drugs.json", "w") as f:
        json.dump({"drugs": drugs}, f, indent=2)
    
    print("Done! JSON files generated in data/")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        # Analyze a single drug from command line argument
        drug_name = sys.argv[1]
        try:
            analysis = analyze_drug(drug_name)
            output_dir = os.path.join(os.path.dirname(__file__), 'data')
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"{drug_name.lower()}.json")
            
            with open(output_file, "w") as f:
                json.dump(analysis, f, indent=2)
            print(f"Analysis complete. Results saved to {output_file}")
            sys.exit(0)
        except Exception as e:
            print(f"Error analyzing drug {drug_name}: {str(e)}", file=sys.stderr)
            sys.exit(1)
    else:
        # Analyze all drugs if no argument provided
        main()
