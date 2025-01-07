import os
import sys
import csv
import re
import time
import requests
import plotly.graph_objects as go
from Bio import Entrez, Medline

# Configure PubMed API
Entrez.email = "your.email@example.com"

def search_pubmed(drug_name):
    """Search PubMed for QT/HR related papers about a drug."""
    query = f'({drug_name}[Title/Abstract]) AND ("QT"[Title/Abstract] OR "heart rate"[Title/Abstract])'
    print(f"\nSearching PubMed with query: {query}")
    
    try:
        print("Performing esearch...")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
        record = Entrez.read(handle)
        handle.close()
        
        print(f"Found {len(record.get('IdList', []))} papers")
        
        papers = []
        if record['IdList']:
            print(f"Paper IDs found: {record['IdList']}")
            print("Fetching paper details...")
            
            handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
            records = list(Medline.parse(handle))
            handle.close()
            
            print(f"Successfully retrieved {len(records)} paper details")
            
            for record in records:
                if 'TI' in record and 'AB' in record:
                    # Get full text if available
                    full_text = get_paper_full_text(record.get('PMID', ''))
                    text_content = full_text if full_text else record.get('AB', '')
                    
                    papers.append({
                        'title': record.get('TI', ''),
                        'text': text_content,
                        'drug_name': drug_name
                    })
                    print(f"\nFound paper: {record.get('TI', '')[:100]}...")
                    print(f"Text length: {len(text_content)} characters")
        return papers
    except Exception as e:
        print(f"PubMed search error: {e}")
        return []

def get_herg_ic50(drug_name):
    """Get hERG IC50 from CSV."""
    try:
        with open('hergic50valueschembl.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['drug_name'].lower() == drug_name.lower():
                    return float(row['herg_ic50'])
    except Exception as e:
        print(f"Error reading hERG IC50: {e}")
    return None

def calculate_concentrations(drug_name, dose_mg):
    """Calculate drug concentrations."""
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/MolecularWeight,XLogP/JSON"
    try:
        resp = requests.get(pubchem_url)
        data = resp.json()['PropertyTable']['Properties'][0]
        mol_weight = float(data['MolecularWeight'])  # Convert to float
        
        dose_mmol = dose_mg / mol_weight
        theoretical_conc = dose_mmol * 1000  # Convert to μM
        bioavailable_conc = theoretical_conc * 0.4  # Assume 40% bioavailability
        herg_ic50 = get_herg_ic50(drug_name) or 0
        
        return {
            'theoretical_conc': theoretical_conc,
            'bioavailable_conc': bioavailable_conc,
            'herg_theoretical_ratio': herg_ic50/theoretical_conc if theoretical_conc else 0,
            'herg_bioavailable_ratio': herg_ic50/bioavailable_conc if bioavailable_conc else 0
        }
    except Exception as e:
        print(f"Error calculating concentrations: {e}")
        return {}

def extract_qt_hr(text):
    """Extract QT and HR values from text."""
    # More comprehensive patterns
    qt_patterns = [
        r"QT[c\s]*[\s:]+(\d{2,3})",  # Basic QT pattern
        r"QT[c\s]*\s*interval[\s:]+(\d{2,3})",  # QT interval
        r"QT[c\s]*\s*prolongation[\s:]*(\d{2,3})",  # QT prolongation
        r"QT[c\s]*\s*of\s+(\d{2,3})",  # QT of X
        r"prolonged\s+QT[c\s]*\s*to\s+(\d{2,3})",  # Prolonged QT to X
        r"QT[c\s]*\s*duration\s+of\s+(\d{2,3})",  # QT duration of X
        r"QT[c\s]*\s*interval\s+was\s+(\d{2,3})",  # QT interval was X
        r"QT[c\s]*\s*interval\s+of\s+(\d{2,3})",  # QT interval of X
        r"QT[c\s]*\s*intervals?\s+(?:ranging\s+)?from\s+\d+\s*(?:to|[-])\s*(\d{2,3})",  # QT interval range (take upper)
        r"QT[c\s]*\s*prolongation\s+(?:up\s+)?to\s+(\d{2,3})"  # QT prolongation up to X
    ]
    
    hr_patterns = [
        r"heart rate[\s:]+(\d{2,3})",  # Basic heart rate
        r"HR[\s:]+(\d{2,3})",  # HR abbreviation
        r"heart rate\s+of\s+(\d{2,3})",  # Heart rate of X
        r"HR\s+of\s+(\d{2,3})",  # HR of X
        r"(\d{2,3})\s*beats\s*(?:per|/)\s*min",  # X beats per minute
        r"(\d{2,3})\s*bpm",  # X bpm
        r"heart rate\s+was\s+(\d{2,3})",  # Heart rate was X
        r"heart rate\s+(?:ranging\s+)?from\s+\d+\s*(?:to|[-])\s*(\d{2,3})",  # Heart rate range (take upper)
        r"bradycardia\s+(?:of|to)\s+(\d{2,3})",  # Bradycardia of/to X
        r"tachycardia\s+(?:of|to)\s+(\d{2,3})"  # Tachycardia of/to X
    ]
    
    # Try all patterns
    qt_value = None
    hr_value = None
    
    text = text.lower()
    print("\nSearching text for QT/HR values:")
    print(text[:500] + "...")  # Print first 500 chars for debugging
    
    # Find QT value
    for pattern in qt_patterns:
        match = re.search(pattern, text)
        if match:
            qt_value = int(match.group(1))
            print(f"Found QT value: {qt_value} (pattern: {pattern})")
            break
    
    # Find HR value
    for pattern in hr_patterns:
        match = re.search(pattern, text)
        if match:
            hr_value = int(match.group(1))
            print(f"Found HR value: {hr_value} (pattern: {pattern})")
            break
    
    if not qt_value and not hr_value:
        print("No QT/HR values found in text")
    
    return {
        'qt': qt_value,
        'hr': hr_value
    }

def get_paper_full_text(pmid):
    """Get full text from PubMed Central if available."""
    try:
        # First check if paper is in PMC
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        record = Entrez.read(handle)
        handle.close()
        
        if record[0]['LinkSetDb']:
            pmc_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
            print(f"Found PMC ID: {pmc_id}")
            
            # Get full text
            handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="text")
            full_text = handle.read()
            handle.close()
            
            return full_text
    except Exception as e:
        print(f"Error getting full text: {e}")
    return None

def create_analysis(drug_name):
    """Create combined analysis with table and nomogram."""
    # Get papers and extract data
    papers = search_pubmed(drug_name)
    data_points = []
    
    for paper in papers:
        values = extract_qt_hr(paper['text'])
        if values['qt'] and values['hr']:
            data_points.append({
                'title': paper['title'],
                'qt': values['qt'],
                'hr': values['hr']
            })
    
    print(f"\nFound {len(data_points)} QT/HR data points from papers")
    
    # Create table
    table = go.Figure(data=[go.Table(
        header=dict(
            values=['Case Report Title', 'QT (ms)', 'Heart Rate (bpm)'],
            font=dict(size=12, color='white'),
            fill_color='royalblue',
            align='left'
        ),
        cells=dict(
            values=[
                [p['title'] for p in data_points],
                [p['qt'] for p in data_points],
                [p['hr'] for p in data_points]
            ],
            font=dict(size=11),
            align='left'
        )
    )])
    
    # Create nomogram
    nomogram = go.Figure()
    if data_points:
        nomogram.add_trace(go.Scatter(
            x=[p['hr'] for p in data_points],
            y=[p['qt'] for p in data_points],
            mode='markers',
            text=[p['title'] for p in data_points],
            hoverinfo='text+x+y'
        ))
    nomogram.update_layout(
        title=f"QT-HR Nomogram for {drug_name}",
        xaxis_title="Heart Rate (bpm)",
        yaxis_title="QT Interval (ms)"
    )
    
    # Save to HTML
    output_file = f"{drug_name}_analysis.html"
    with open(output_file, 'w') as f:
        f.write("<html><body>")
        f.write(table.to_html(full_html=False))
        f.write(nomogram.to_html(full_html=False))
        f.write("</body></html>")
    print(f"Saved combined analysis to: {os.path.abspath(output_file)}")

def main():
    if len(sys.argv) < 2:
        print("Please provide a drug name")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    print(f"\nAnalyzing {drug_name}...")
    
    # Calculate concentrations for standard doses
    test_doses = [25, 50]
    for dose in test_doses:
        print(f"\nCalculating concentrations for {dose}mg dose:")
        results = calculate_concentrations(drug_name, dose)
        if results:
            for key, value in results.items():
                print(f"{key}: {value:.3f}")
    
    # Create analysis
    create_analysis(drug_name)
    print(f"\nAnalysis complete! Open {drug_name}_analysis.html to view results.")

if __name__ == "__main__":
    main()
