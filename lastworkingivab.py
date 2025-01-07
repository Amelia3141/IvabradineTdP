import os
import sys
import csv
import re
import time
import requests
import plotly.graph_objects as go
from Bio import Entrez, Medline
from typing import Dict, List, Optional, Tuple, Union
import xml.etree.ElementTree as ET

# Configure PubMed API
Entrez.email = "ghhercock@gmail.com"

def search_pubmed(drug_name: str) -> List[Dict[str, str]]:
    """Search PubMed for QT/HR related papers about a drug."""
    query = f'({drug_name}[Title/Abstract]) AND ("QT"[Title/Abstract] OR "heart rate"[Title/Abstract] OR "torsade"[Title/Abstract])'
    print(f"\nSearching PubMed with query: {query}")
    
    try:
        # First search
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle)
        handle.close()
        
        papers = []
        if record['IdList']:
            print(f"Found {len(record['IdList'])} papers")
            print(f"Paper IDs found: {record['IdList']}")
            
            # Get paper details
            handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
            records = list(Medline.parse(handle))
            handle.close()
            
            print(f"Successfully retrieved {len(records)} paper details\n")
            
            for record in records:
                if 'TI' in record:
                    # Try to get full text
                    full_text = get_paper_full_text(record)
                    text_length = len(record.get('AB', '')) if record.get('AB') else 0
                    
                    print(f"\nFound paper: {record.get('TI', '')[:100]}...")
                    print(f"Text length: {text_length} characters")
                    
                    papers.append({
                        'title': record.get('TI', ''),
                        'abstract': record.get('AB', ''),
                        'full_text': full_text,
                        'pmid': record.get('PMID', ''),
                        'doi': record.get('DOI', ''),
                        'year': record.get('DP', '')[:4] if record.get('DP') else ''
                    })
        
        return papers
    except Exception as e:
        print(f"PubMed search error: {e}")
        return []

def get_paper_full_text(paper_info):
    """Try to get full text from multiple sources."""
    pmid = paper_info.get('PMID', '')
    doi = paper_info.get('DOI', '')
    title = paper_info.get('TI', '')
    
    # 1. Try PubMed Central first
    try:
        if pmid:
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
            record = Entrez.read(handle)
            handle.close()
            
            if record[0].get('LinkSetDb', []):
                pmc_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
                print(f"Found PMC ID: {pmc_id}")
                
                try:
                    # Try to get XML format first for better parsing
                    handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml", retmode="xml")
                    xml_content = handle.read()
                    handle.close()
                    
                    # Extract text from XML, focusing on methods, results, and discussion sections
                    text_content = []
                    xml_root = ET.fromstring(xml_content)
                    
                    # Extract abstract
                    abstract = xml_root.find(".//abstract")
                    if abstract is not None:
                        text_content.append(ET.tostring(abstract, encoding='unicode', method='text'))
                    
                    # Extract main sections
                    for section in xml_root.findall(".//sec"):
                        title_elem = section.find("title")
                        if title_elem is not None:
                            title_text = title_elem.text.lower()
                            if any(word in title_text for word in ['method', 'result', 'discussion', 'case', 'observation']):
                                text_content.append(ET.tostring(section, encoding='unicode', method='text'))
                    
                    if text_content:
                        return ' '.join(text_content)
                    
                    # If XML parsing fails, fall back to text format
                    handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="text")
                    full_text = handle.read()
                    handle.close()
                    return full_text
                    
                except Exception as e:
                    print(f"Error getting full text: {e}")
    except Exception as e:
        print(f"PMC retrieval error: {e}")
    
    # 2. Try bioRxiv/medRxiv
    try:
        if doi:
            biorxiv_url = f"https://api.biorxiv.org/details/biorxiv/{doi}"
            medrxiv_url = f"https://api.biorxiv.org/details/medrxiv/{doi}"
            
            for url in [biorxiv_url, medrxiv_url]:
                response = requests.get(url)
                if response.ok:
                    data = response.json()
                    if data.get('collection', []):
                        pdf_url = data['collection'][0].get('pdf_url')
                        if pdf_url:
                            print(f"Found preprint PDF: {pdf_url}")
                            return f"PDF available at: {pdf_url}"
    except Exception as e:
        print(f"Preprint retrieval error: {e}")
    
    return None

def get_herg_ic50(drug_name: str) -> Optional[float]:
    """Get hERG IC50 from CSV."""
    try:
        with open('hergic50valueschembl.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['Drug'].lower() == drug_name.lower():
                    return float(row['hERG IC50 (μM)'])
    except Exception as e:
        print(f"Error reading hERG IC50: {e}")
    return None

def calculate_concentrations(dose_mg, molecular_weight=427.5, volume_of_distribution=252, bioavailability=0.4):
    """Calculate theoretical and bioavailable concentrations."""
    try:
        # Convert dose from mg to μmol
        dose_umol = (dose_mg / molecular_weight) * 1000
        
        # Calculate theoretical concentration (μM)
        theoretical_conc = dose_umol / volume_of_distribution
        
        # Calculate bioavailable concentration (μM)
        bioavailable_conc = theoretical_conc * bioavailability
        
        return {
            'theoretical_conc': round(theoretical_conc, 3),
            'bioavailable_conc': round(bioavailable_conc, 3)
        }
    except Exception as e:
        print(f"Error calculating concentrations: {e}")
        return {
            'theoretical_conc': 0,
            'bioavailable_conc': 0
        }

def fetch_pubchem_data(drug_name: str) -> Optional[Dict[str, str]]:
    """
    Fetch molecular weight and other properties from PubChem.
    Returns None if data not found or error occurs.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/MolecularWeight,XLogP/JSON"
    
    try:
        print(f"Fetching PubChem data for {drug_name}...")
        response = requests.get(url)
        response.raise_for_status()
        
        data = response.json()['PropertyTable']['Properties'][0]
        print(f"Found PubChem CID: {data.get('CID')}")
        print(f"Molecular Weight: {data.get('MolecularWeight')} g/mol")
        print(f"XLogP: {data.get('XLogP')}")
        
        return data
        
    except Exception as e:
        print(f"Error fetching PubChem data: {e}")
        return None

def calculate_drug_concentrations(drug_name: str, dose_mg: float) -> Dict[str, Union[float, str]]:
    """Calculate drug concentrations and hERG ratios using PubChem data and known hERG values."""
    if not dose_mg:
        return {
            'theoretical_conc': 'N/A',
            'bioavailable_conc': 'N/A',
            'herg_theoretical_ratio': 'N/A',
            'herg_bioavailable_ratio': 'N/A'
        }
    
    print(f"\nCalculating concentrations for {dose_mg}mg dose of {drug_name}...")
    
    # Fetch drug data from PubChem
    pubchem_data = fetch_pubchem_data(drug_name)
    if not pubchem_data:
        return {
            'theoretical_conc': 'Unknown',
            'bioavailable_conc': 'Unknown',
            'herg_theoretical_ratio': 'Unknown',
            'herg_bioavailable_ratio': 'Unknown'
        }
    
    # Get hERG IC50 value from CSV
    herg_ic50 = get_herg_ic50(drug_name)
    
    # Use PubChem data
    params = {
        'molecular_weight': float(pubchem_data['MolecularWeight']),  # Convert to float
        'herg_ic50': herg_ic50,  # Will be None if not found
        'bioavailability': 0.5  # Default to 50% if not available
    }
    
    # Estimate bioavailability based on XLogP (lipophilicity)
    # Higher XLogP generally correlates with better oral bioavailability
    xlogp = pubchem_data.get('XLogP')
    if xlogp is not None:
        if xlogp < 0:
            params['bioavailability'] = 0.2  # Poor absorption
        elif xlogp < 2:
            params['bioavailability'] = 0.4  # Moderate absorption
        elif xlogp < 4:
            params['bioavailability'] = 0.6  # Good absorption
        else:
            params['bioavailability'] = 0.8  # Excellent absorption
        print(f"Estimated bioavailability: {params['bioavailability']*100:.0f}% based on XLogP of {xlogp}")
    
    # Calculate concentrations
    dose_mmol = dose_mg / params['molecular_weight'] * 1000  # Convert to mmol
    theoretical_conc_um = dose_mmol * 1000  # Convert to μM
    bioavailable_conc_um = theoretical_conc_um * params['bioavailability']
    
    print(f"Dose in mmol: {dose_mmol:.3f}")
    print(f"Theoretical concentration: {theoretical_conc_um:.2f} μM")
    print(f"Bioavailable concentration: {bioavailable_conc_um:.2f} μM")
    
    # Calculate ratios if hERG IC50 is available
    if params['herg_ic50'] is not None:
        herg_theoretical_ratio = params['herg_ic50'] / theoretical_conc_um
        herg_bioavailable_ratio = params['herg_ic50'] / bioavailable_conc_um
        herg_theoretical_ratio_str = f"{herg_theoretical_ratio:.3f}"
        herg_bioavailable_ratio_str = f"{herg_bioavailable_ratio:.3f}"
        print(f"hERG IC50/Theoretical Ratio: {herg_theoretical_ratio_str}")
        print(f"hERG IC50/Bioavailable Ratio: {herg_bioavailable_ratio_str}")
    else:
        herg_theoretical_ratio_str = 'N/A'
        herg_bioavailable_ratio_str = 'N/A'
        print("No hERG IC50 value available for ratio calculations")
    
    return {
        'theoretical_conc': f"{theoretical_conc_um:.2f}",
        'bioavailable_conc': f"{bioavailable_conc_um:.2f}",
        'herg_theoretical_ratio': herg_theoretical_ratio_str,
        'herg_bioavailable_ratio': herg_bioavailable_ratio_str
    }

def calculate_naranjo_score(text: str) -> int:
    """Calculate Naranjo adverse drug reaction probability score."""
    score = 0
    
    # Previous conclusive reports
    if re.search(r'previous.*report|similar.*case|known.*reaction', text, re.I):
        score += 1
    
    # Adverse event after drug
    if re.search(r'after.*administration|following.*drug|post.*treatment', text, re.I):
        score += 2
    
    # Improvement after discontinuation
    if re.search(r'improve.*discontinu|resolve.*stop|better.*withdraw', text, re.I):
        score += 1
    
    # Recurrence on readministration
    if re.search(r'recur.*readminist|return.*restart|reappear.*rechallenge', text, re.I):
        score += 2
    
    # Alternative causes
    if not re.search(r'other.*cause|alternative.*explanation|different.*reason', text, re.I):
        score += 2
    
    # Reaction with placebo
    if re.search(r'placebo', text, re.I):
        score -= 1
    
    # Drug levels in blood
    if re.search(r'blood.*level|plasma.*concentration|serum.*drug', text, re.I):
        score += 1
    
    # Dose-response relationship
    if re.search(r'dose.*response|concentration.*dependent|dose.*related', text, re.I):
        score += 1
    
    # Previous exposure
    if re.search(r'previous.*exposure|prior.*treatment|history.*medication', text, re.I):
        score += 1
    
    # Objective evidence
    if re.search(r'ECG|EKG|QT|QTc|rhythm|arrhythmia|monitor', text, re.I):
        score += 1
    
    return score

def calculate_tisdale_score(text: str) -> int:
    """Calculate Tisdale risk score for QT prolongation."""
    score = 0
    
    # Age ≥68 years
    if re.search(r'(?:age|year).{0,10}(?:6[8-9]|[7-9][0-9]|1[0-9][0-9])', text, re.I):
        score += 1
    
    # Female sex
    if re.search(r'female|woman|women', text, re.I):
        score += 1
    
    # Loop diuretic
    if re.search(r'furosemide|bumetanide|torsemide|ethacrynic acid', text, re.I):
        score += 1
    
    # Serum K+ ≤3.5
    if re.search(r'potassium.{0,20}(?:1\.[0-9]|2\.[0-9]|3\.[0-4])|K\+.{0,10}(?:1\.[0-9]|2\.[0-9]|3\.[0-4])', text, re.I):
        score += 2
    
    # Admission QTc ≥450 ms
    if re.search(r'QTc.{0,20}(?:4[5-9][0-9]|[5-9][0-9][0-9])', text, re.I):
        score += 2
    
    # Acute MI
    if re.search(r'myocardial infarction|MI|heart attack', text, re.I):
        score += 2
    
    # Sepsis
    if re.search(r'sepsis|septic|infection', text, re.I):
        score += 3
    
    # Heart failure
    if re.search(r'heart failure|CHF|cardiac failure|reduced ejection', text, re.I):
        score += 3
    
    # One QTc-prolonging drug
    if re.search(r'antiarrhythmic|antipsychotic|antibiotic|antidepressant', text, re.I):
        score += 3
    
    # ≥2 QTc-prolonging drugs
    if len(re.findall(r'antiarrhythmic|antipsychotic|antibiotic|antidepressant', text, re.I)) >= 2:
        score += 3
    
    return score

def assess_who_umc(text: str) -> str:
    """Assess WHO-UMC causality categories."""
    # Certain
    if (re.search(r'rechallenge.*positive|readministration.*recur', text, re.I) and
        re.search(r'dechallenge.*positive|discontinuation.*improve', text, re.I) and
        re.search(r'plausible.*time|temporal.*relation', text, re.I)):
        return "Certain"
    
    # Probable
    elif (re.search(r'dechallenge.*positive|discontinuation.*improve', text, re.I) and
          re.search(r'plausible.*time|temporal.*relation', text, re.I) and
          not re.search(r'alternative.*cause|other.*explanation', text, re.I)):
        return "Probable"
    
    # Possible
    elif (re.search(r'plausible.*time|temporal.*relation', text, re.I) and
          re.search(r'alternative.*cause|other.*explanation', text, re.I)):
        return "Possible"
    
    # Unlikely
    elif re.search(r'implausible.*time|unlikely.*relation|improbable', text, re.I):
        return "Unlikely"
    
    # Conditional
    elif re.search(r'more.*data|additional.*information|further.*investigation', text, re.I):
        return "Conditional"
    
    # Unassessable
    else:
        return "Unassessable"

def save_to_csv(data: List[Dict[str, str]], filepath: str) -> None:
    """Save list of dictionaries to CSV file."""
    if not data:
        return
    
    fieldnames = data[0].keys()
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def process_excel_data(excel_path: str) -> List[Dict[str, float]]:
    """Process data from CSV file instead of Excel."""
    data_points = []
    
    with open(excel_path, 'r') as file:
        reader = csv.DictReader(file)
        print("Raw data from CSV:")
        for row in reader:
            try:
                qt = float(row['QT'])
                hr = float(row['HR'])
                data_points.append({'hr': hr, 'qt': qt})
                print(f"Successfully added point: HR={hr}, QT={qt}")
            except (ValueError, KeyError) as e:
                print(f"Error processing row: {row}, Error: {str(e)}")
    
    print(f"\nTotal points processed: {len(data_points)}")
    return data_points

def extract_values_from_text(text: str) -> List[Dict[str, float]]:
    """Extract QT intervals and heart rates from text using NLP."""
    # Load English language model
    nlp = spacy.load("en_core_web_sm")
    
    # Process the text
    doc = nlp(text)
    
    # Patterns for QT and HR values
    qt_patterns = [
        r"(?:QT|QTc)\s*(?:interval|duration)?\s*(?:was|of|:)?\s*(\d+)\s*(?:ms|msec|milliseconds)",
        r"(?:uncorrected\s+)?QT\s+interval\s+(?:of|was|:)?\s*(\d+)\s*(?:ms|msec|milliseconds)"
    ]
    
    hr_patterns = [
        r"(?:heart\s+rate|HR)\s*(?:was|of|:)?\s*(\d+)\s*(?:bpm|beats per minute)",
        r"(?:heart\s+rate|HR)\s*(?:ranging|ranged)\s*(?:from)?\s*(\d+)\s*(?:to|-)?\s*\d+\s*(?:bpm|beats per minute)"
    ]
    
    qt_values = []
    hr_values = []
    
    # Find all matches in the text
    for pattern in qt_patterns:
        matches = re.finditer(pattern, text.lower())
        for match in matches:
            try:
                value = int(match.group(1))
                if 200 <= value <= 1000:  # Reasonable QT range in ms
                    qt_values.append(value)
            except ValueError:
                continue

    for pattern in hr_patterns:
        matches = re.finditer(pattern, text.lower())
        for match in matches:
            try:
                value = int(match.group(1))
                if 20 <= value <= 300:  # Reasonable HR range
                    hr_values.append(value)
            except ValueError:
                continue
    
    # Pair QT and HR values if we have equal numbers of each
    data_points = []
    for qt, hr in zip(qt_values, hr_values):
        data_points.append({'qt': qt, 'hr': hr})
    
    return data_points

def process_pdf(pdf_path: str) -> List[Dict[str, float]]:
    """Process a PDF file and extract QT and HR values."""
    data_points = []
    
    with open(pdf_path, 'rb') as file:
        reader = PyPDF2.PdfReader(file)
        
        # Process each page
        for page in reader.pages:
            text = page.extract_text()
            if text:
                page_data = extract_values_from_text(text)
                data_points.extend(page_data)
    
    # Save to CSV
    if data_points:
        csv_path = pdf_path.replace('.pdf', '.csv')
        save_to_csv(data_points, csv_path)
        print(f"\nSaved extracted data to: {csv_path}")
        print(f"\nTotal points found: {len(data_points)}")
    
    return data_points

def process_multiple_pdfs(pdf_directory: str) -> List[Dict[str, float]]:
    """Process all PDFs in a directory."""
    all_data_points = []
    pdf_files = list(Path(pdf_directory).glob('*.pdf'))
    
    print(f"Found {len(pdf_files)} PDF files in directory")
    
    for pdf_path in pdf_files:
        print(f"\nProcessing: {pdf_path.name}")
        data_points = process_pdf(pdf_path)
        all_data_points.extend(data_points)
    
    return all_data_points

def extract_case_details(text: str) -> Dict[str, Union[str, int, float, None]]:
    """Extract detailed case information using NLP patterns."""
    # Load English language model
    nlp = spacy.load("en_core_web_sm")
    doc = nlp(text.lower())
    
    # Initialize case details dictionary
    case = {
        'Case Report Title': '',
        'Age': None,
        'Sex': None,
        'Oral Dose': None,
        'theoretical max concentration': None,
        '40% bioavailability': None,
        'Theoretical HERG IC50': None,
        '40% Plasma concentration': None,
        'Uncorrected QT': None,
        'QTB': None,
        'QTF': None,
        'Heart Rate': None,
        'Torsades de Pointes?': 'N',
        'Blood Pressure': None,
        'Medical History': '',
        'Medication History': '',
        'Course Treatment': '',
        'Naranjo Score': None,
        'WHO-UMC': '',
        'Tisdale Score': None
    }
    
    # Extract age
    age_pattern = r"(?:age[d]?|year[s]?[-\s]old)\D*(\d+)"
    age_matches = re.findall(age_pattern, text.lower())
    if age_matches:
        case['Age'] = int(age_matches[0])
    
    # Extract sex
    if re.search(r"\b(fe)?male\b", text.lower()):
        case['Sex'] = 'Female' if 'female' in text.lower() else 'Male'
    
    # Extract dose
    dose_pattern = r"(\d+(?:\.\d+)?)\s*(?:mg|milligrams?)(?:\s*oral)?"
    dose_matches = re.findall(dose_pattern, text.lower())
    if dose_matches:
        case['Oral Dose'] = float(dose_matches[0])
    
    # Extract QT intervals
    qt_pattern = r"(?:uncorrected\s+)?qt\s*(?:interval|duration)?\s*(?:was|of|:)?\s*(\d+)\s*(?:ms|msec|milliseconds)"
    qt_matches = re.findall(qt_pattern, text.lower())
    if qt_matches:
        case['Uncorrected QT'] = int(qt_matches[0])
    
    # Extract heart rate
    hr_pattern = r"(?:heart\s+rate|hr)\s*(?:was|of|:)?\s*(\d+)\s*(?:bpm|beats per minute)"
    hr_matches = re.findall(hr_pattern, text.lower())
    if hr_matches:
        case['Heart Rate'] = int(hr_matches[0])
    
    # Extract blood pressure
    bp_pattern = r"(?:blood\s+pressure|bp)\s*(?:was|of|:)?\s*(\d+/\d+)"
    bp_matches = re.findall(bp_pattern, text.lower())
    if bp_matches:
        case['Blood Pressure'] = bp_matches[0]
    
    # Check for TdP
    if re.search(r"torsades\s+de\s+pointes|tdp", text.lower()):
        case['Torsades de Pointes?'] = 'Y'
    
    # Extract medical history (look for section or keywords)
    med_history_pattern = r"(?:medical history|past medical history|pmh)[:\.](.*?)(?:\n|\.)"
    med_history_matches = re.findall(med_history_pattern, text.lower())
    if med_history_matches:
        case['Medical History'] = med_history_matches[0].strip()
    
    # Extract medication history
    med_pattern = r"(?:medication[s]?|drugs|treatment)[:\.](.*?)(?:\n|\.)"
    med_matches = re.findall(med_pattern, text.lower())
    if med_matches:
        case['Medication History'] = med_matches[0].strip()
    
    # Extract Naranjo and Tisdale scores if present
    naranjo_pattern = r"naranjo\s*(?:score|scale)?\s*(?:was|of|:)?\s*(\d+)"
    naranjo_matches = re.findall(naranjo_pattern, text.lower())
    if naranjo_matches:
        case['Naranjo Score'] = int(naranjo_matches[0])
    
    tisdale_pattern = r"tisdale\s*(?:score|scale)?\s*(?:was|of|:)?\s*(\d+)"
    tisdale_matches = re.findall(tisdale_pattern, text.lower())
    if tisdale_matches:
        case['Tisdale Score'] = int(tisdale_matches[0])
    
    return case

def process_pdf_for_cases(pdf_path: str) -> List[Dict[str, Union[str, int, float, None]]]:
    """Process a PDF file and extract case report details."""
    cases = []
    
    with open(pdf_path, 'rb') as file:
        reader = PyPDF2.PdfReader(file)
        
        # Process each page
        full_text = ""
        for page in reader.pages:
            full_text += page.extract_text() + "\n"
        
        # Extract case details
        case = extract_case_details(full_text)
        case['Case Report Title'] = os.path.basename(pdf_path).replace('.pdf', '')
        cases.append(case)
    
    # Save to CSV
    csv_path = pdf_path.replace('.pdf', '_cases.csv')
    save_to_csv(cases, csv_path)
    print(f"\nSaved case analysis to: {csv_path}")
    
    return cases

def process_multiple_pdfs_for_cases(pdf_directory: str) -> List[Dict[str, Union[str, int, float, None]]]:
    """Process all PDFs in a directory and compile case details."""
    all_cases = []
    pdf_files = list(Path(pdf_directory).glob('*.pdf'))
    
    print(f"Processing {len(pdf_files)} PDF files for case details")
    
    for pdf_path in pdf_files:
        print(f"\nProcessing: {pdf_path.name}")
        cases = process_pdf_for_cases(pdf_path)
        all_cases.extend(cases)
    
    return all_cases

def get_reference_lines() -> List[Dict[str, Union[str, List[Dict[str, float]], str]]]:
    """Get reference lines for QT nomogram."""
    # Solid line coordinates
    nomogram_line_solid = [
        { 'hr': 10, 'qt': 480 },
        { 'hr': 60, 'qt': 480 },
        { 'hr': 110, 'qt': 400 }
    ]
    
    # Dashed line coordinates - extend the last point to ensure even dashes
    nomogram_line_dashed = [
        { 'hr': 110, 'qt': 400 },
        { 'hr': 145, 'qt': 330 },
        { 'hr': 150, 'qt': 322 }  # Added extra point to extend the line
    ]
    
    return [
        {'data': nomogram_line_solid, 'name': 'High Risk', 'color': 'red', 'dash': 'solid'},
        {'data': nomogram_line_dashed, 'name': 'Normal Upper Limit', 'color': 'blue', 'dash': 'dash'}
    ]

def plot_qt_nomogram(data_points: List[Dict[str, float]]) -> go.Figure:
    """Plot QT nomogram using Plotly."""
    # Create reference lines
    ref_lines = get_reference_lines()
    
    # Create figure
    fig = go.Figure()
    
    # Add reference lines
    for line in ref_lines:
        hr_values = [point['hr'] for point in line['data']]
        qt_values = [point['qt'] for point in line['data']]
        fig.add_trace(go.Scatter(
            x=hr_values,
            y=qt_values,
            mode='lines',
            name=line['name'],
            line=dict(color=line['color'], dash=line['dash']),
            showlegend=True
        ))
    
    # Add data points
    if data_points:
        hr_values = [point['hr'] for point in data_points]
        qt_values = [point['qt'] for point in data_points]
        fig.add_trace(go.Scatter(
            x=hr_values,
            y=qt_values,
            mode='markers',
            name='Data Points',
            marker=dict(color='red')
        ))
    
    # Update layout
    fig.update_layout(
        title='QT Nomogram',
        xaxis_title='Heart Rate (bpm)',
        yaxis_title='QT Interval (ms)',
        showlegend=True
    )
    
    return fig

def extract_dose(text: str) -> Optional[float]:
    """Extract drug dose from text."""
    # Look for dose patterns like "5 mg", "5mg", "5-10mg", etc.
    dose_pattern = r'(\d+(?:\.\d+)?)\s*(?:mg|milligrams?)(?:\s*oral)?'
    match = re.search(dose_pattern, text.lower())
    if match:
        return float(match.group(1))
    return None

def plot_papers_table(papers: List[Dict[str, str]]) -> go.Figure:
    """Create a table showing paper details and QT/HR data."""
    # Prepare table headers
    headers = ['Case Report Title', 'Age', 'Sex', 'Oral Dose (mg)', 
              'theoretical max concentration (μM)', '40% bioavailability Plasma concentration μM',
              'Theoretical HERG IC50 / Concentration μM', '40% Plasma concentration/HERG IC50',
              'Uncorrected QT (ms)', 'QTc', 'QTB', 'QTF', 'Heart Rate (bpm)',
              'Torsades de Pointes?', 'Blood Pressure (mmHg)', 'Medical History',
              'Medication History', 'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC']
    
    # Extract and analyze data
    analyzed_data = analyze_papers_for_qt_hr(papers)
    
    # Prepare table data
    data = []
    seen_papers = set()
    
    for point in analyzed_data:
        if point['title'] not in seen_papers:
            row = [
                point['title'],  # Case Report Title
                '',  # Age
                '',  # Sex
                '',  # Oral Dose
                '',  # Theoretical max concentration
                '',  # 40% bioavailability
                '',  # Theoretical HERG IC50
                '',  # 40% Plasma concentration/HERG IC50
                point['qt'],  # Uncorrected QT
                '',  # QTc
                '',  # QTB
                '',  # QTF
                point['hr'],  # Heart Rate
                '',  # Torsades
                '',  # Blood Pressure
                '',  # Medical History
                '',  # Medication History
                '',  # Course Treatment
                point['naranjo_score'],  # Naranjo Score
                point['who_umc'],  # WHO-UMC
                point['tisdale_score']  # Tisdale Score
            ]
            data.append(row)
            seen_papers.add(point['title'])
    
    # Create table with updated styling
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=headers,
            font=dict(size=11, color='white'),
            fill_color='royalblue',
            align=['left'] * len(headers),
            height=40
        ),
        cells=dict(
            values=list(zip(*data)) if data else [[] for _ in headers],
            align=['left'] * len(headers),
            font=dict(size=10),
            height=30,
            fill=dict(color=['rgb(245, 245, 245)', 'white'])
        )
    )])
    
    # Update layout with better sizing
    fig.update_layout(
        title=dict(
            text=f"Literature Analysis Results for {papers[0]['drug_name'] if papers else 'Unknown Drug'}",
            x=0.5,
            font=dict(size=20)
        ),
        width=1500,  # Increased width to accommodate more columns
        height=max(len(data) * 40 + 100, 400),  # Dynamic height based on rows
        margin=dict(t=50, l=20, r=20, b=20)
    )
    
    return fig

def get_paper_abstract(pmid: str) -> Optional[str]:
    """Get paper abstract from PubMed."""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        return abstract
    except Exception as e:
        print(f"Error fetching abstract for {pmid}: {e}")
        return None

def plot_combined_analysis(papers: List[Dict[str, str]], drug_name: str, qt_data: List[Dict[str, float]]) -> None:
    """Create a combined figure with QT nomogram and analysis table."""
    # Create subplots with more space for the table
    fig = sp.make_subplots(
        rows=2, cols=1,
        row_heights=[0.6, 0.4],  # Give more height to the table
        specs=[[{"type": "scatter"}],
               [{"type": "table"}]],
        vertical_spacing=0.1,
        subplot_titles=(f"QT Nomogram for {drug_name}", f"Literature Analysis for {drug_name}")
    )

    # Add reference lines
    ref_lines = get_reference_lines()
    for line in ref_lines:
        hr_values = [point['hr'] for point in line['data']]
        qt_values = [point['qt'] for point in line['data']]
        fig.add_trace(
            go.Scatter(
                x=hr_values,
                y=qt_values,
                mode='lines',
                name=line['name'],
                line=dict(color=line['color'], dash=line['dash']),
                showlegend=True
            ),
            row=1, col=1
        )

    # Add data points
    if qt_data:
        hr_values = [point['hr'] for point in qt_data]
        qt_values = [point['qt'] for point in qt_data]
        fig.add_trace(
            go.Scatter(
                x=hr_values,
                y=qt_values,
                mode='markers',
                name='Data Points',
                marker=dict(color='red')
            ),
            row=1, col=1
        )

    # Update nomogram layout
    fig.update_xaxes(title_text='Heart Rate (bpm)', range=[0, 300], row=1, col=1)
    fig.update_yaxes(title_text='QT Interval (ms)', range=[200, 800], row=1, col=1)

    # Create table data with more columns
    table_data = {
        'Title': [paper['Title'] for paper in papers],
        'DOI': [paper['DOI'] for paper in papers],
        'Year': [paper.get('Year', 'N/A') for paper in papers],
        'QT Values': [paper.get('QT_values', 'N/A') for paper in papers],
        'HR Values': [paper.get('HR_values', 'N/A') for paper in papers],
        'Naranjo': [paper.get('Naranjo_Score', 'N/A') for paper in papers],
        'Tisdale': [paper.get('Tisdale_Score', 'N/A') for paper in papers],
        'WHO-UMC': [paper.get('WHO_UMC', 'N/A') for paper in papers]
    }

    # Add table with improved formatting
    fig.add_trace(
        go.Table(
            header=dict(
                values=list(table_data.keys()),
                fill_color='paleturquoise',
                align='left',
                font=dict(size=12, color='black'),
                height=40
            ),
            cells=dict(
                values=list(table_data.values()),
                fill_color='lavender',
                align='left',
                font=dict(size=11),
                height=30,
                line_color='darkslategray',
                line_width=1
            )
        ),
        row=2, col=1
    )

    # Update layout
    fig.update_layout(
        height=1500,  # Increased height to accommodate table
        showlegend=True,
        title_text=f"Combined Analysis for {drug_name}",
        title_x=0.5,
        paper_bgcolor='white',
        plot_bgcolor='aliceblue'
    )

    fig.show()

def process_paper(paper: Dict[str, str]) -> Dict[str, Union[str, int, float, None]]:
    """Process a single paper and extract all relevant information."""
    paper_info = {
        'Title': paper.get('title', 'N/A'),
        'DOI': paper.get('doi', 'N/A'),
        'Authors': paper.get('authors', 'N/A'),
        'Age': 'N/A',
        'Sex': 'N/A',
        'QT_values': 'N/A',
        'HR_values': 'N/A',
        'BP_values': 'N/A',
        'Medical_History': 'N/A',
        'Naranjo_Score': 0,
        'Tisdale_Score': 0,
        'WHO_UMC': 'N/A'
    }
    
    text = paper.get('abstract', '') + ' ' + paper.get('title', '')
    if not text:
        return paper_info
    
    # Extract age using improved patterns
    age_patterns = [
        r'(\d+)[-\s]year[-\s]old',
        r'age[d]?\s+(\d+)',
        r'(\d+)\s*(?:yo|y\.o\.|years\s+old)',
        r'(\d+)\s*(?:year|yr)[-\s]old',
        r'aged?\s*(?:of)?\s*(\d+)\s*(?:year|yr)',
    ]
    
    for pattern in age_patterns:
        match = re.search(pattern, text.lower())
        if match:
            try:
                age = int(match.group(1))
                if 0 <= age <= 120:  # Reasonable age range
                    paper_info['Age'] = age
                    break
            except ValueError:
                continue
    
    # Extract sex
    sex_patterns = [
        r'\b(male|female)\b',
        r'\b(man|woman)\b',
        r'\b(boy|girl)\b',
        r'\b(gentleman|lady)\b'
    ]
    
    for pattern in sex_patterns:
        match = re.search(pattern, text.lower())
        if match:
            sex = match.group(1)
            if sex in ['male', 'man', 'boy', 'gentleman']:
                paper_info['Sex'] = 'Male'
            elif sex in ['female', 'woman', 'girl', 'lady']:
                paper_info['Sex'] = 'Female'
            break
    
    # Extract clinical values
    clinical_values = extract_qt_hr_values(text)
    
    # Format values for CSV
    if clinical_values['qt_values']:
        paper_info['QT_values'] = ', '.join(map(str, clinical_values['qt_values']))
    if clinical_values['hr_values']:
        paper_info['HR_values'] = ', '.join(map(str, clinical_values['hr_values']))
    if clinical_values['bp_values']:
        bp_str = []
        for bp in clinical_values['bp_values']:
            if isinstance(bp, tuple) and bp[1] is not None:
                bp_str.append(f"{bp[0]}/{bp[1]}")
            else:
                bp_str.append(str(bp[0]))
        paper_info['BP_values'] = ', '.join(bp_str)
    
    # Extract medical history
    medical_history = extract_medical_history(text)
    if medical_history:
        paper_info['Medical_History'] = medical_history
    
    # Calculate scores
    paper_info['Naranjo_Score'] = calculate_naranjo_score(text)
    paper_info['Tisdale_Score'] = calculate_tisdale_score(text)
    paper_info['WHO_UMC'] = assess_who_umc(text)
    
    return paper_info

def save_to_csv(papers: List[Dict[str, str]], filepath: str) -> None:
    """Save paper analysis to CSV file with all extracted information."""
    if not papers:
        return
    
    fieldnames = [
        'Title', 'Year', 'DOI', 'Authors', 'Age', 'Sex', 'QT_values', 'HR_values', 'BP_values',
        'Medical_History', 'Naranjo_Score', 'Tisdale_Score', 'WHO_UMC', 'TdP_Present'
    ]
    
    processed_papers = []
    for paper in papers:
        paper_info = process_paper(paper)
        
        # Add TdP presence
        text = paper.get('abstract', '') + ' ' + paper.get('title', '')
        paper_info['TdP_Present'] = 'Yes' if any(term in text.lower() for term in ['torsade', 'tdp', 'torsades']) else 'No'
        
        # Add year if available
        paper_info['Year'] = paper.get('year', 'N/A')
        
        processed_papers.append(paper_info)
    
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(processed_papers)

def extract_qt_hr_values(text: str) -> Dict[str, List[Union[int, Tuple[int, str, str]]]]:
    """
    Extract QT, QTc, and HR values from text with context.
    Returns a dictionary containing:
    - qt_values: List of tuples (value, method, context)
    - qtc_values: List of tuples (value, formula, context)
    - hr_values: List of tuples (value, method, context)
    """
    values = {
        'qt_values': [],
        'qtc_values': [],
        'hr_values': [],
        'bp_values': []
    }
    
    # QT interval patterns with improved context capture
    qt_patterns = [
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:of|was|is|=|:|measured|prolonged|increased|decreased)\s*(?:to|at|as)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:duration|measurement|value)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:baseline\s+)?QT\s+interval\s+(?:range|varied|ranging|from)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:prolongation|increase|decrease)\s*(?:to|of|by)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})'
    ]
    
    # QTc patterns with formula detection
    qtc_patterns = [
        r'(?P<pre>[^.]{0,100})QTc\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:of|was|is|=|:|measured)\s*(?:to|at|as)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QTc\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:interval|duration|measurement|value)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:corrected|heart[\s-]rate[\s-]corrected)\s+QT\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QTc\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:prolongation|increase|decrease)\s*(?:to|of|by)?\s*(?P<value>\d{2,3})\s*(?:±\s*\d+)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})'
    ]
    
    # Heart rate patterns with improved context
    hr_patterns = [
        r'(?P<pre>[^.]{0,100})(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:of|was|is|=|:|measured)\s*(?:to|at|as)?\s*(?P<value>\d{1,3})\s*(?:±\s*\d+)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:increased|decreased|changed)\s*(?:to|at)?\s*(?P<value>\d{1,3})\s*(?:±\s*\d+)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:baseline\s+)?(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:range|varied|ranging|from)?\s*(?P<value>\d{1,3})\s*(?:±\s*\d+)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:mean|average)\s+(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{1,3})\s*(?:±\s*\d+)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})'
    ]
    
    def extract_with_context(patterns, text, value_range, value_type):
        results = []
        for pattern in patterns:
            matches = re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE | re.DOTALL)
            for match in matches:
                try:
                    value = int(match.group('value'))
                    min_val, max_val = value_range
                    
                    if min_val <= value <= max_val:
                        pre_context = match.group('pre').strip()
                        post_context = match.group('post').strip()
                        
                        # Clean up the context
                        context = f"{pre_context} [...] {post_context}"
                        context = re.sub(r'\s+', ' ', context).strip()
                        context = context.replace('[...] [...]', '[...]')
                        
                        # Remove any incomplete sentences
                        context = re.sub(r'^[a-z].*?\.\.\.', '', context)
                        context = re.sub(r'\.\.\.[a-z].*$', '', context)
                        
                        if value_type == 'qtc' and 'formula' in match.groupdict():
                            formula = match.group('formula') or 'unspecified'
                            results.append((value, formula, context))
                        else:
                            method = 'measured' if any(word in pre_context.lower() for word in ['measured', 'recorded', 'observed', 'mean', 'average']) else 'reported'
                            results.append((value, method, context))
                except (ValueError, IndexError):
                    continue
        return results

    # Extract values with context
    values['qt_values'] = extract_with_context(qt_patterns, text, (200, 1000), 'qt')
    values['qtc_values'] = extract_with_context(qtc_patterns, text, (200, 1000), 'qtc')
    values['hr_values'] = extract_with_context(hr_patterns, text, (20, 300), 'hr')
    
    # Remove duplicates while preserving context
    values['qt_values'] = list(dict.fromkeys(values['qt_values']))
    values['qtc_values'] = list(dict.fromkeys(values['qtc_values']))
    values['hr_values'] = list(dict.fromkeys(values['hr_values']))
    
    return values

def extract_dose(text: str) -> Optional[float]:
    """Extract drug dose from text."""
    dose_patterns = [
        r"(\d+)\s*mg(?:\s*of\s*\w+)?",
        r"dose\s*(?:of\s*)?(\d+)\s*mg",
        r"(\d+)\s*milligrams?"
    ]
    
    for pattern in dose_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            try:
                return int(match.group(1))
            except ValueError:
                continue
    return None

def extract_medical_history(text: str) -> str:
    """Extract medical history from text."""
    conditions = []
    
    # Common conditions and their variations
    condition_patterns = [
        r'(?:history of|diagnosed with|known|chronic|suffers? from)\s+([\w\s-]+(?:failure|disease|disorder|syndrome|attack|arrhythmia|hypertension|diabetes|cancer|infection))',
        r'(?:underlying|pre-existing)\s+([\w\s-]+(?:condition|disease|disorder))',
        r'(?:chronic|acute)\s+([\w\s-]+(?:failure|disease|disorder))',
        r'(\w+\s+(?:failure|disease|disorder))',
        r'((?:dilated|ischemic|alcoholic)\s+cardiomyopathy)',
        r'(heart\s+failure(?:\s+with\s+\w+\s+ejection\s+fraction)?)',
        r'(\w+\s+kidney\s+disease)',
        r'(coronary\s+artery\s+disease)',
        r'(atrial\s+fibrillation)',
        r'((?:type\s+[12]|gestational)\s+diabetes)',
        r'(hypertension)',
        r'(myocardial\s+infarction)',
        r'(cardiac\s+arrest)',
        r'((?:tachy|brady)(?:cardia|arrhythmia))',
        r'((?:prolonged|long)\s+QT[a-z]*\s+(?:syndrome)?)',
        r'(asthma|copd|bronchitis)',
        r'(respiratory\s+(?:failure|distress|infection))',
        r'(pneumonia)',
        r'(seizure|epilepsy)',
        r'(stroke|tia)',
        r'(parkinson)',
        r'(depression)',
        r'(anxiety)',
        r'(bipolar\s+disorder)',
        r'(schizophrenia)',
        r'(diabetes(?:\s+mellitus)?)',
        r'(renal\s+(?:failure|insufficiency))',
        r'(liver\s+(?:failure|disease|cirrhosis))',
        r'(cancer|malignancy)',
        r'(sepsis|infection)',
        r'(covid-19|coronavirus)',
        r'((?:low|high)\s+(?:potassium|magnesium|calcium))',
        r'(electrolyte\s+(?:imbalance|abnormality))',
    ]
    
    for pattern in condition_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            condition = match.group(1).strip()
            if condition and len(condition) > 3:  # Avoid short matches
                conditions.append(condition)
    
    # Remove duplicates while preserving order
    conditions = list(dict.fromkeys(conditions))
    return '; '.join(conditions) if conditions else 'N/A'

def extract_medications(text: str) -> str:
    """Extract medication history from text."""
    medications = []
    
    # Common drug classes and specific drugs
    drug_patterns = [
        r'(?:taking|on|receiving|treated with|therapy with)\s+([\w\s-]+(?:blocker|inhibitor|antagonist|agonist|sartan|pril|statin))',
        r'(?:prescribed|administered|given|received)\s+([\w\s-]+(?:mg|mcg|g|ml|units))',
        r'(\w+)\s+(?:\d+\s*(?:mg|mcg|g|ml|units))',
        r'((?:beta|alpha)\s*-?\s*blockers?)',
        r'(ace\s*inhibitors?)',
        r'(calcium\s*channel\s*blockers?)',
        r'(diuretics?)',
        r'(anti-?\s*arrhythmics?)',
        r'(anti-?\s*coagulants?)',
        r'(anti-?\s*depressants?)',
        r'(anti-?\s*psychotics?)',
        r'(anti-?\s*biotics?)',
        r'(anti-?\s*virals?)',
        r'\b(warfarin|heparin|aspirin|clopidogrel)',
        r'\b(metoprolol|atenolol|carvedilol|bisoprolol)',
        r'\b(lisinopril|enalapril|ramipril)',
        r'\b(amlodipine|diltiazem|verapamil)',
        r'\b(furosemide|bumetanide|spironolactone)',
        r'\b(amiodarone|sotalol|flecainide)',
        r'\b(digoxin)',
        r'\b(azithromycin|clarithromycin|erythromycin)',
        r'\b(ciprofloxacin|levofloxacin|moxifloxacin)',
        r'\b(penicillin|amoxicillin|ampicillin)',
        r'\b(ceftriaxone|cefuroxime|cephalexin)',
        r'\b(fluoxetine|sertraline|paroxetine|citalopram|escitalopram)',
        r'\b(quetiapine|olanzapine|risperidone|haloperidol)',
        r'\b(lithium|valproate|carbamazepine)',
        r'\b(metformin|insulin)',
        r'\b(hydroxychloroquine|chloroquine)',
        r'\b(ivermectin|remdesivir)',
    ]
    
    for pattern in drug_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            medication = match.group(1).strip()
            if medication and len(medication) > 2:  # Avoid short matches
                medications.append(medication)
    
    # Remove duplicates while preserving order
    medications = list(dict.fromkeys(medications))
    return '; '.join(medications) if medications else 'N/A'

def extract_treatment_course(text: str) -> str:
    """Extract treatment course from text."""
    treatments = []
    
    # Treatment patterns
    treatment_patterns = [
        r'(?:treated|managed)\s+with\s+([\w\s-]+)',
        r'(?:administered|given|received)\s+([\w\s-]+)',
        r'(?:therapy|treatment)\s+with\s+([\w\s-]+)',
        r'(defibrillation|cardioversion|pacing)',
        r'(magnesium\s+(?:infusion|supplementation))',
        r'(potassium\s+(?:infusion|supplementation))',
        r'(activated\s+charcoal)',
        r'(intubation|ventilation)',
        r'(cpr|cardiopulmonary\s+resuscitation)',
        r'(hemodialysis|dialysis)',
        r'(ecg\s+monitoring)',
        r'(telemetry)',
        r'(discontinuation|withdrawal|stopped)\s+of\s+([\w\s-]+)',
        r'(dose\s+reduction|decreased|increased)\s+of\s+([\w\s-]+)',
        r'(monitoring|observation)',
        r'(iv\s+fluids?)',
        r'(oxygen\s+therapy)',
        r'(mechanical\s+ventilation)',
        r'(vasopressors?)',
        r'(inotropes?)',
        r'(antibiotics?)',
        r'(antiarrhythmics?)',
    ]
    
    for pattern in treatment_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            if len(match.groups()) > 1:
                treatment = f"{match.group(1)} {match.group(2)}".strip()
            else:
                treatment = match.group(1).strip()
            if treatment and len(treatment) > 3:
                treatments.append(treatment)
    
    # Remove duplicates while preserving order
    treatments = list(dict.fromkeys(treatments))
    return '; '.join(treatments) if treatments else 'N/A'

def extract_qt_hr_with_context(text: str) -> List[Dict[str, Union[int, str]]]:
    """Extract QT and HR values with their surrounding context."""
    qt_data = []
    
    # Split text into sentences for better context
    sentences = text.split('.')
    for sentence in sentences:
        qt, qt_match = extract_qtc(sentence)
        hr, hr_match = extract_hr(sentence)
        
        if qt and hr:
            # Basic physiological validation
            if 200 <= qt <= 800 and 20 <= hr <= 300:
                qt_data.append({
                    'qt': qt,
                    'hr': hr,
                    'context': sentence.strip(),
                    'qt_match': qt_match,
                    'hr_match': hr_match
                })
    return qt_data

def analyze_papers_for_qt_hr(papers: List[Dict[str, str]]) -> List[Dict[str, Union[int, str]]]:
    """Analyze papers to extract QT/HR data with context."""
    analyzed_data = []
    
    for paper in papers:
        text = paper['title'] + ' ' + paper['abstract']
        qt_hr_points = extract_qt_hr_with_context(text)
        
        if qt_hr_points:
            for point in qt_hr_points:
                analyzed_data.append({
                    'title': paper['title'],
                    'doi': paper.get('doi', 'N/A'),
                    'qt': point['qt'],
                    'hr': point['hr'],
                    'context': point['context'],
                    'qt_match': point['qt_match'],
                    'hr_match': point['hr_match'],
                    'naranjo_score': paper['naranjo_score'],
                    'tisdale_score': paper['tisdale_score'],
                    'who_umc': paper['who_umc']
                })
                print(f"\nFound QT/HR data point in paper:")
                print(f"Title: {paper['title']}")
                print(f"QT: {point['qt']} ms (matched: {point['qt_match']})")
                print(f"HR: {point['hr']} bpm (matched: {point['hr_match']})")
                print(f"Context: {point['context']}")
    
    return analyzed_data

def create_qt_nomogram(drug_name: str) -> go.Figure:
    """Create QT nomogram using data points from analyzed papers."""
    # Create figure with white background
    fig = go.Figure(layout=dict(
        plot_bgcolor='white',
        paper_bgcolor='white'
    ))
    
    # Get papers and analyze for QT/HR data
    papers = search_pubmed(drug_name)
    analyzed_data = analyze_papers_for_qt_hr(papers)
    
    print(f"\nFound {len(analyzed_data)} QT/HR data points from papers")
    
    # Add data points from papers
    if analyzed_data:
        hr_values = [point['hr'] for point in analyzed_data]
        qt_values = [point['qt'] for point in analyzed_data]
        hover_texts = [
            f"Paper: {point['title']}<br>"
            f"QT: {point['qt']} ms<br>"
            f"HR: {point['hr']} bpm<br>"
            f"Context: {point['context']}"
            for point in analyzed_data
        ]
        
        fig.add_trace(go.Scatter(
            x=hr_values,
            y=qt_values,
            mode='markers',
            name='Data Points',
            text=hover_texts,
            hoverinfo='text'
        ))
    
    # Add nomogram lines
    hr_range = list(range(30, 121, 1))
    
    # Upper line (high risk)
    upper_qt = [440 + 0.675 * (60 - hr) for hr in hr_range]
    fig.add_trace(go.Scatter(
        x=hr_range,
        y=upper_qt,
        mode='lines',
        name='High Risk',
        line=dict(color='red', dash='solid', width=2)
    ))
    
    # Lower line (normal)
    lower_qt = [400 + 0.675 * (60 - hr) for hr in hr_range]
    fig.add_trace(go.Scatter(
        x=hr_range,
        y=lower_qt,
        mode='lines',
        name='Normal Limit',
        line=dict(color='blue', dash='dash', width=2)
    ))
    
    # Update layout with grid and better styling
    fig.update_layout(
        title=dict(
            text=f"QT Nomogram for {drug_name}",
            x=0.5,
            font=dict(size=20)
        ),
        xaxis=dict(
            title="Heart Rate (bpm)",
            gridcolor='lightgray',
            showgrid=True,
            zeroline=False,
            range=[0, 300]
        ),
        yaxis=dict(
            title="QT Interval (ms)",
            gridcolor='lightgray',
            showgrid=True,
            zeroline=False,
            range=[200, 800]
        ),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        ),
        width=1000,
        height=800,
        hovermode='closest'
    )
    
    return fig

def plot_papers_table(papers: List[Dict[str, str]]) -> go.Figure:
    """Create a table showing paper details and QT/HR data."""
    # Prepare table headers
    headers = [
        'Case Report Title', 'Age', 'Sex', 'Oral Dose (mg)', 
        'theoretical max concentration (μM)', '40% bioavailability Plasma concentration μM',
        'Theoretical HERG IC50 / Concentration μM', '40% Plasma concentration/HERG IC50',
        'Uncorrected QT (ms)', 'QTc', 'QTB', 'QTF', 'Heart Rate (bpm)',
        'Torsades de Pointes?', 'Blood Pressure (mmHg)', 'Medical History',
        'Medication History', 'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC'
    ]
    
    # Extract and analyze data
    analyzed_data = analyze_papers_for_qt_hr(papers)
    
    # Prepare table data
    data = []
    seen_papers = set()
    
    for point in analyzed_data:
        if point['title'] not in seen_papers:
            row = [
                point['title'],  # Case Report Title
                '',  # Age
                '',  # Sex
                '',  # Oral Dose
                '',  # Theoretical max concentration
                '',  # 40% bioavailability
                '',  # Theoretical HERG IC50
                '',  # 40% Plasma concentration/HERG IC50
                point['qt'],  # Uncorrected QT
                '',  # QTc
                '',  # QTB
                '',  # QTF
                point['hr'],  # Heart Rate
                '',  # Torsades
                '',  # Blood Pressure
                '',  # Medical History
                '',  # Medication History
                '',  # Course Treatment
                point['naranjo_score'],  # Naranjo Score
                point['who_umc'],  # WHO-UMC
                point['tisdale_score']  # Tisdale Score
            ]
            data.append(row)
            seen_papers.add(point['title'])
    
    # Create table with updated styling
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=headers,
            font=dict(size=11, color='white'),
            fill_color='royalblue',
            align=['left'] * len(headers),
            height=40
        ),
        cells=dict(
            values=list(zip(*data)) if data else [[] for _ in headers],
            align=['left'] * len(headers),
            font=dict(size=10),
            height=30,
            fill=dict(color=['rgb(245, 245, 245)', 'white'])
        )
    )])
    
    # Update layout with better sizing
    fig.update_layout(
        title=dict(
            text=f"Literature Analysis Results for {papers[0]['drug_name'] if papers else 'Unknown Drug'}",
            x=0.5,
            font=dict(size=20)
        ),
        width=1500,  # Increased width to accommodate more columns
        height=max(len(data) * 40 + 100, 400),  # Dynamic height based on rows
        margin=dict(t=50, l=20, r=20, b=20)
    )
    
    return fig

def create_analysis(drug_name: str) -> None:
    """Create combined analysis with table and nomogram."""
    # Get papers and extract data
    papers = search_pubmed(drug_name)
    data_points = []
    
    for paper in papers:
        values = extract_qt_hr(paper['abstract'])
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
            name='Data Points',
            text=[p['title'] for p in data_points],
            hoverinfo='text'
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

def extract_qt_hr(text: str) -> Dict[str, Optional[int]]:
    """Extract QT and HR values from text."""
    qt_pattern = r"QT[c\s]*[\s:]+(\d{2,3})"
    hr_pattern = r"heart rate[\s:]+(\d{2,3})"
    
    qt = re.search(qt_pattern, text.lower())
    hr = re.search(hr_pattern, text.lower())
    
    return {
        'qt': int(qt.group(1)) if qt else None,
        'hr': int(hr.group(1)) if hr else None
    }

def main(drug_name: Optional[str] = None) -> None:
    """Main function to run the TdP risk analysis."""
    if drug_name is None:
        if len(sys.argv) < 2:
            print("Please provide a drug name as argument")
            sys.exit(1)
        drug_name = sys.argv[1]
    
    print(f"\nAnalyzing {drug_name.title()}...")
    print("--------------------------------\n")
    
    # Test concentration calculations with a standard dose
    test_doses = [25, 50]  # Common doses in mg
    for dose in test_doses:
        print(f"\nCalculating concentrations for {dose}mg dose:")
        results = calculate_concentrations(dose)
        print(f"Theoretical concentration: {results['theoretical_conc']} μM")
        print(f"Bioavailable concentration: {results['bioavailable_conc']} μM")
    
    # Create output directory
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Create combined analysis
    print("\nCreating combined analysis...")
    create_analysis(drug_name)
    print(f"\nAnalysis complete! Open {drug_name}_analysis.html to view results.")

if __name__ == "__main__":
    main()
