import os
import plotly.graph_objects as go
import plotly.subplots as sp
import numpy as np
import PyPDF2
import spacy
import re
import csv
from pathlib import Path
from Bio import Entrez
import requests
import time
from urllib.parse import urlencode
import json
import argparse
from io import StringIO
from Bio import Medline
import logging

# Set your email for PubMed API
Entrez.email = "ghhercock@gmail.com"  # Replace with your email

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    nlp = spacy.load('en_core_web_sm')
except OSError:
    import subprocess
    subprocess.run(['python', '-m', 'spacy', 'download', 'en_core_web_sm'])
    nlp = spacy.load('en_core_web_sm')

def calculate_naranjo_score(text):
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

def calculate_tisdale_score(text):
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

def assess_who_umc(text):
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

def save_to_csv(data, filepath):
    """Save list of dictionaries to CSV file."""
    if not data:
        return
    
    fieldnames = data[0].keys()
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def process_excel_data(excel_path):
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

def extract_values_from_text(text):
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

def process_pdf(pdf_path):
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

def process_multiple_pdfs(pdf_directory):
    """Process all PDFs in a directory."""
    all_data_points = []
    pdf_files = list(Path(pdf_directory).glob('*.pdf'))
    
    print(f"Found {len(pdf_files)} PDF files in directory")
    
    for pdf_path in pdf_files:
        print(f"\nProcessing: {pdf_path.name}")
        data_points = process_pdf(pdf_path)
        all_data_points.extend(data_points)
    
    return all_data_points

def extract_case_details(text):
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

def process_pdf_for_cases(pdf_path):
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

def process_multiple_pdfs_for_cases(pdf_directory):
    """Process all PDFs in a directory and compile case details."""
    all_cases = []
    pdf_files = list(Path(pdf_directory).glob('*.pdf'))
    
    print(f"Processing {len(pdf_files)} PDF files for case details")
    
    for pdf_path in pdf_files:
        print(f"\nProcessing: {pdf_path.name}")
        cases = process_pdf_for_cases(pdf_path)
        all_cases.extend(cases)
    
    return all_cases

def get_reference_lines():
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

def plot_qt_nomogram(data_points):
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

def calculate_drug_concentrations(drug_name, dose_mg):
    """Calculate drug concentrations and hERG ratios."""
    if not dose_mg:
        return {
            'theoretical_conc': 'N/A',
            'bioavailable_conc': 'N/A',
            'herg_theoretical_ratio': 'N/A',
            'herg_bioavailable_ratio': 'N/A'
        }
    
    # Drug-specific parameters
    drug_params = {
        'Ivabradine': {
            'molecular_weight': 468.593,  # g/mol
            'herg_ic50': 3.0,  # μM
            'bioavailability': 0.40  # 40%
        },
        'Citalopram': {
            'molecular_weight': 324.39,  # g/mol
            'herg_ic50': 1.5,  # μM
            'bioavailability': 0.80  # 80%
        }
    }
    
    params = drug_params.get(drug_name)
    if not params:
        return {
            'theoretical_conc': 'Unknown',
            'bioavailable_conc': 'Unknown',
            'herg_theoretical_ratio': 'Unknown',
            'herg_bioavailable_ratio': 'Unknown'
        }
    
    # Calculate concentrations
    dose_mmol = dose_mg / params['molecular_weight'] * 1000  # Convert to mmol
    theoretical_conc_um = dose_mmol * 1000  # Convert to μM
    bioavailable_conc_um = theoretical_conc_um * params['bioavailability']
    
    # Calculate ratios
    herg_theoretical_ratio = params['herg_ic50'] / theoretical_conc_um
    herg_bioavailable_ratio = params['herg_ic50'] / bioavailable_conc_um
    
    return {
        'theoretical_conc': f"{theoretical_conc_um:.2f}",
        'bioavailable_conc': f"{bioavailable_conc_um:.2f}",
        'herg_theoretical_ratio': f"{herg_theoretical_ratio:.3f}",
        'herg_bioavailable_ratio': f"{herg_bioavailable_ratio:.3f}"
    }

def extract_dose(text):
    """Extract drug dose from text."""
    # Look for dose patterns like "5 mg", "5mg", "5-10mg", etc.
    dose_pattern = r'(\d+(?:\.\d+)?)\s*(?:mg|milligrams?)(?:\s*oral)?'
    match = re.search(dose_pattern, text.lower())
    if match:
        return float(match.group(1))
    return None

def plot_papers_table(papers, drug_name):
    """Create a table showing paper details and causality scores."""
    # Extract relevant data
    data = []
    columns = ['Title', 'Year', 'Dose (mg)', 'Theoretical Conc (μM)', 
              'Bioavailable Conc (μM)', 'hERG IC50/Theoretical', 
              'hERG IC50/Bioavailable', 'QTc (ms)', 'HR (bpm)',
              'Naranjo', 'Tisdale', 'WHO-UMC']
    
    for paper in papers:
        # Extract dose from title and abstract
        text = paper['title'] + ' ' + paper['abstract']
        
        # Calculate drug concentrations
        conc_data = calculate_drug_concentrations(drug_name, extract_dose(text))
        
        # Extract QTc and HR
        qtc_match = re.search(r'(?:QTc|QT).*?(\d+)\s*(?:ms|msec)', text)
        hr_match = re.search(r'(?:HR|heart rate).*?(\d+)\s*(?:bpm|beats)', text)
        
        row = [
            paper['title'][:50] + '...' if len(paper['title']) > 50 else paper['title'],
            paper['year'],
            str(extract_dose(text)) if extract_dose(text) else 'N/A',
            conc_data['theoretical_conc'],
            conc_data['bioavailable_conc'],
            conc_data['herg_theoretical_ratio'],
            conc_data['herg_bioavailable_ratio'],
            qtc_match.group(1) if qtc_match else 'N/A',
            hr_match.group(1) if hr_match else 'N/A',
            str(paper['naranjo_score']),
            str(paper['tisdale_score']),
            paper['who_umc']
        ]
        data.append(row)
    
    # Create figure and axis
    fig = go.Figure(data=[go.Table(header=dict(values=columns),
                                   cells=dict(values=data))])
    
    # Update layout
    fig.update_layout(title=f'{drug_name} Case Reports Analysis')
    
    # Save figure
    fig.write_html(f'{drug_name.lower()}_analysis_table.html')
    fig.write_pdf(f'{drug_name.lower()}_analysis_table.pdf')
    return fig

def get_paper_abstract(pmid):
    """Get paper abstract from PubMed."""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        return abstract
    except Exception as e:
        print(f"Error fetching abstract for {pmid}: {e}")
        return None

def search_pubmed(mesh_terms, title_abs_terms):
    """Search PubMed for papers matching the terms."""
    papers = []
    try:
        # Search using MeSH terms
        handle = Entrez.esearch(db="pubmed", term=mesh_terms)
        record = Entrez.read(handle)
        mesh_ids = record["IdList"]
        
        # Search using title/abstract terms
        handle = Entrez.esearch(db="pubmed", term=title_abs_terms)
        record = Entrez.read(handle)
        title_abs_ids = record["IdList"]
        
        # Combine and deduplicate IDs
        pmids = list(set(mesh_ids + title_abs_ids))
        
        # Fetch details for each paper
        for pmid in pmids:
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                records = list(Medline.parse(handle))
                if records:
                    record = records[0]
                    if 'TI' in record:  # Title exists
                        # Clean up DOI - it might be in different fields
                        doi = None
                        if 'LID' in record:
                            doi_match = re.search(r'(10\.\d{4,}/\S+)', record['LID'])
                            if doi_match:
                                doi = doi_match.group(1)
                        if not doi and 'AID' in record:
                            for id in record['AID']:
                                if '[doi]' in id:
                                    doi = id.replace(' [doi]', '')
                                    break
                        
                        paper = {
                            'title': record.get('TI', ''),
                            'abstract': record.get('AB', ''),
                            'authors': ', '.join(record.get('AU', [])),
                            'journal': record.get('JT', ''),
                            'year': record.get('DP', '').split()[0] if record.get('DP') else '',
                            'doi': doi or 'N/A',
                            'naranjo_score': calculate_naranjo_score(record.get('AB', '')),
                            'tisdale_score': calculate_tisdale_score(record.get('AB', '')),
                            'who_umc': assess_who_umc(record.get('AB', ''))
                        }
                        papers.append(paper)
                        print(f"\nFound paper: {paper['title']}")
                        print(f"DOI: {paper['doi']}")
                        print("Causality Assessment:")
                        print(f"- Naranjo Score: {paper['naranjo_score']}")
                        print(f"- Tisdale Score: {paper['tisdale_score']}")
                        print(f"- WHO-UMC: {paper['who_umc']}\n")
            except Exception as e:
                print(f"Error processing paper {pmid}: {e}")
                continue
                
    except Exception as e:
        print(f"Error searching PubMed: {e}")
    
    return papers

def save_to_csv(papers, output_path):
    """Save paper analysis to CSV file."""
    headers = ['Title', 'Year', 'DOI', 'Authors', 'Age', 'Sex', 'Dose (mg)', 
              'QT (ms)', 'QTc (ms)', 'QTB (ms)', 'QTF (ms)', 'HR (bpm)', 
              'TdP?', 'Medical History', 'Medication History',
              'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC']
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        
        for paper in papers:
            text = paper['title'] + ' ' + paper['abstract']
            
            # Extract age and sex
            age_match = re.search(r'(\d+)[\s-]*(?:year|yr)', text, re.IGNORECASE)
            sex_match = re.search(r'(male|female|man|woman|boy|girl)', text, re.IGNORECASE)
            
            row = {
                'Title': paper['title'],
                'Year': paper['year'],
                'DOI': paper.get('doi', 'N/A'),
                'Authors': paper['authors'],
                'Age': age_match.group(1) if age_match else 'N/A',
                'Sex': sex_match.group(1).title() if sex_match else 'N/A',
                'Dose (mg)': extract_dose(text) or 'N/A',
                'QT (ms)': extract_qt_hr_values(text)['qt_values'][0] if extract_qt_hr_values(text)['qt_values'] else 'N/A',
                'QTc (ms)': extract_qt_hr_values(text)['qt_values'][1] if len(extract_qt_hr_values(text)['qt_values']) > 1 else 'N/A',
                'QTB (ms)': extract_qt_hr_values(text)['qt_values'][2] if len(extract_qt_hr_values(text)['qt_values']) > 2 else 'N/A',
                'QTF (ms)': extract_qt_hr_values(text)['qt_values'][3] if len(extract_qt_hr_values(text)['qt_values']) > 3 else 'N/A',
                'HR (bpm)': extract_qt_hr_values(text)['hr_values'][0] if extract_qt_hr_values(text)['hr_values'] else 'N/A',
                'TdP?': 'Yes' if 'torsade' in text.lower() else 'N/A',
                'Medical History': extract_medical_history(text),
                'Medication History': extract_medications(text),
                'Course of Treatment': extract_treatment_course(text),
                'Outcome': extract_outcome(text),
                'Naranjo': str(paper['naranjo_score']),
                'Tisdale': str(paper['tisdale_score']),
                'WHO-UMC': paper['who_umc']
            }
            writer.writerow(row)

def search_pubmed_papers(drug_name):
    """Search PubMed for case reports of TdP with the specified drug."""
    # MESH term query - replace drug names
    mesh_query = f'((hERG) OR (QT) OR (QTc) OR (torsad*)) AND ({drug_name})'
    
    # Additional title/abstract terms
    title_abs_query = f'({drug_name}[Title/Abstract]) AND (("torsades de pointes"[Title/Abstract]) OR (tdp[Title/Abstract]) OR (torsades[Title/Abstract]))'
    
    # Combine queries with case report filter
    query = f'({mesh_query} OR {title_abs_query}) AND (case reports[Publication Type])'
    
    print(f"\nSearching PubMed using:")
    print(f"MESH terms: {mesh_query}")
    print(f"Title/Abstract terms: {title_abs_query}")
    
    # Search PubMed
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        print("No papers found matching the criteria.")
        return []
    
    print(f"\nFound {len(record['IdList'])} matching papers")
    
    # Fetch details for each paper
    handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="xml", retmode="xml")
    papers = Entrez.read(handle)
    handle.close()
    
    paper_details = []
    for paper in papers['PubmedArticle']:
        try:
            article = paper['MedlineCitation']['Article']
            
            # Get DOI if available
            doi = None
            if 'ELocationID' in article:
                for id in article['ELocationID']:
                    if id.attributes['EIdType'] == 'doi':
                        doi = str(id)
                        break
            
            # Get full text if available
            full_text = ""
            if 'Abstract' in article:
                if isinstance(article['Abstract']['AbstractText'], list):
                    full_text = ' '.join([str(text) for text in article['Abstract']['AbstractText']])
                else:
                    full_text = str(article['Abstract']['AbstractText'])
            
            paper_info = {
                'pmid': paper['MedlineCitation']['PMID'],
                'title': article['ArticleTitle'],
                'authors': ', '.join([author['LastName'] + ' ' + author['ForeName'] for author in article['AuthorList']]),
                'journal': article['Journal']['Title'],
                'year': article['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'abstract': full_text,
                'doi': doi
            }
            
            # Calculate causality scores
            combined_text = f"{paper_info['title']} {paper_info['abstract']}"
            paper_info['naranjo_score'] = calculate_naranjo_score(combined_text)
            paper_info['tisdale_score'] = calculate_tisdale_score(combined_text)
            paper_info['who_umc'] = assess_who_umc(combined_text)
            
            # Print detailed information
            print(f"\nFound paper: {paper_info['title']}")
            print(f"DOI: {paper_info['doi'] if paper_info['doi'] else 'N/A'}")
            print(f"Causality Assessment:")
            print(f"- Naranjo Score: {paper_info['naranjo_score']}")
            print(f"- Tisdale Score: {paper_info['tisdale_score']}")
            print(f"- WHO-UMC: {paper_info['who_umc']}")
            
            paper_details.append(paper_info)
            
        except Exception as e:
            print(f"Error processing paper: {str(e)}")
    
    return paper_details

def plot_combined_analysis(papers, drug_name, qt_data):
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
                line=dict(
                    color=line['color'],
                    dash=line['dash']
                )
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
                marker=dict(
                    size=10,
                    color='red',
                    symbol='circle'
                )
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

    return fig

def get_paper_abstract(pmid):
    """Get paper abstract from PubMed."""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        return abstract
    except Exception as e:
        print(f"Error fetching abstract for {pmid}: {e}")
        return None

def search_pubmed(mesh_terms, title_abs_terms):
    """Search PubMed for papers matching the terms."""
    papers = []
    try:
        # Search using MeSH terms
        handle = Entrez.esearch(db="pubmed", term=mesh_terms)
        record = Entrez.read(handle)
        mesh_ids = record["IdList"]
        
        # Search using title/abstract terms
        handle = Entrez.esearch(db="pubmed", term=title_abs_terms)
        record = Entrez.read(handle)
        title_abs_ids = record["IdList"]
        
        # Combine and deduplicate IDs
        pmids = list(set(mesh_ids + title_abs_ids))
        
        # Fetch details for each paper
        for pmid in pmids:
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                records = list(Medline.parse(handle))
                if records:
                    record = records[0]
                    if 'TI' in record:  # Title exists
                        # Clean up DOI - it might be in different fields
                        doi = None
                        if 'LID' in record:
                            doi_match = re.search(r'(10\.\d{4,}/\S+)', record['LID'])
                            if doi_match:
                                doi = doi_match.group(1)
                        if not doi and 'AID' in record:
                            for id in record['AID']:
                                if '[doi]' in id:
                                    doi = id.replace(' [doi]', '')
                                    break
                        
                        paper = {
                            'title': record.get('TI', ''),
                            'abstract': record.get('AB', ''),
                            'authors': ', '.join(record.get('AU', [])),
                            'journal': record.get('JT', ''),
                            'year': record.get('DP', '').split()[0] if record.get('DP') else '',
                            'doi': doi or 'N/A',
                            'naranjo_score': calculate_naranjo_score(record.get('AB', '')),
                            'tisdale_score': calculate_tisdale_score(record.get('AB', '')),
                            'who_umc': assess_who_umc(record.get('AB', ''))
                        }
                        papers.append(paper)
                        print(f"\nFound paper: {paper['title']}")
                        print(f"DOI: {paper['doi']}")
                        print("Causality Assessment:")
                        print(f"- Naranjo Score: {paper['naranjo_score']}")
                        print(f"- Tisdale Score: {paper['tisdale_score']}")
                        print(f"- WHO-UMC: {paper['who_umc']}\n")
            except Exception as e:
                print(f"Error processing paper {pmid}: {e}")
                continue
                
    except Exception as e:
        print(f"Error searching PubMed: {e}")
    
    return papers

def save_to_csv(papers, output_path):
    """Save paper analysis to CSV file."""
    headers = ['Title', 'Year', 'DOI', 'Authors', 'Age', 'Sex', 'Dose (mg)', 
              'QT (ms)', 'QTc (ms)', 'QTB (ms)', 'QTF (ms)', 'HR (bpm)', 
              'TdP?', 'Medical History', 'Medication History',
              'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC']
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        
        for paper in papers:
            text = paper['title'] + ' ' + paper['abstract']
            
            # Extract age and sex
            age_match = re.search(r'(\d+)[\s-]*(?:year|yr)', text, re.IGNORECASE)
            sex_match = re.search(r'(male|female|man|woman|boy|girl)', text, re.IGNORECASE)
            
            row = {
                'Title': paper['title'],
                'Year': paper['year'],
                'DOI': paper.get('doi', 'N/A'),
                'Authors': paper['authors'],
                'Age': age_match.group(1) if age_match else 'N/A',
                'Sex': sex_match.group(1).title() if sex_match else 'N/A',
                'Dose (mg)': extract_dose(text) or 'N/A',
                'QT (ms)': extract_qt_hr_values(text)['qt_values'][0] if extract_qt_hr_values(text)['qt_values'] else 'N/A',
                'QTc (ms)': extract_qt_hr_values(text)['qt_values'][1] if len(extract_qt_hr_values(text)['qt_values']) > 1 else 'N/A',
                'QTB (ms)': extract_qt_hr_values(text)['qt_values'][2] if len(extract_qt_hr_values(text)['qt_values']) > 2 else 'N/A',
                'QTF (ms)': extract_qt_hr_values(text)['qt_values'][3] if len(extract_qt_hr_values(text)['qt_values']) > 3 else 'N/A',
                'HR (bpm)': extract_qt_hr_values(text)['hr_values'][0] if extract_qt_hr_values(text)['hr_values'] else 'N/A',
                'TdP?': 'Yes' if 'torsade' in text.lower() else 'N/A',
                'Medical History': extract_medical_history(text),
                'Medication History': extract_medications(text),
                'Course of Treatment': extract_treatment_course(text),
                'Outcome': extract_outcome(text),
                'Naranjo': str(paper['naranjo_score']),
                'Tisdale': str(paper['tisdale_score']),
                'WHO-UMC': paper['who_umc']
            }
            writer.writerow(row)

def search_pubmed_papers(drug_name):
    """Search PubMed for case reports of TdP with the specified drug."""
    # MESH term query - replace drug names
    mesh_query = f'((hERG) OR (QT) OR (QTc) OR (torsad*)) AND ({drug_name})'
    
    # Additional title/abstract terms
    title_abs_query = f'({drug_name}[Title/Abstract]) AND (("torsades de pointes"[Title/Abstract]) OR (tdp[Title/Abstract]) OR (torsades[Title/Abstract]))'
    
    # Combine queries with case report filter
    query = f'({mesh_query} OR {title_abs_query}) AND (case reports[Publication Type])'
    
    print(f"\nSearching PubMed using:")
    print(f"MESH terms: {mesh_query}")
    print(f"Title/Abstract terms: {title_abs_query}")
    
    # Search PubMed
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        print("No papers found matching the criteria.")
        return []
    
    print(f"\nFound {len(record['IdList'])} matching papers")
    
    # Fetch details for each paper
    handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="xml", retmode="xml")
    papers = Entrez.read(handle)
    handle.close()
    
    paper_details = []
    for paper in papers['PubmedArticle']:
        try:
            article = paper['MedlineCitation']['Article']
            
            # Get DOI if available
            doi = None
            if 'ELocationID' in article:
                for id in article['ELocationID']:
                    if id.attributes['EIdType'] == 'doi':
                        doi = str(id)
                        break
            
            # Get full text if available
            full_text = ""
            if 'Abstract' in article:
                if isinstance(article['Abstract']['AbstractText'], list):
                    full_text = ' '.join([str(text) for text in article['Abstract']['AbstractText']])
                else:
                    full_text = str(article['Abstract']['AbstractText'])
            
            paper_info = {
                'pmid': paper['MedlineCitation']['PMID'],
                'title': article['ArticleTitle'],
                'authors': ', '.join([author['LastName'] + ' ' + author['ForeName'] for author in article['AuthorList']]),
                'journal': article['Journal']['Title'],
                'year': article['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'abstract': full_text,
                'doi': doi
            }
            
            # Calculate causality scores
            combined_text = f"{paper_info['title']} {paper_info['abstract']}"
            paper_info['naranjo_score'] = calculate_naranjo_score(combined_text)
            paper_info['tisdale_score'] = calculate_tisdale_score(combined_text)
            paper_info['who_umc'] = assess_who_umc(combined_text)
            
            # Print detailed information
            print(f"\nFound paper: {paper_info['title']}")
            print(f"DOI: {paper_info['doi'] if paper_info['doi'] else 'N/A'}")
            print(f"Causality Assessment:")
            print(f"- Naranjo Score: {paper_info['naranjo_score']}")
            print(f"- Tisdale Score: {paper_info['tisdale_score']}")
            print(f"- WHO-UMC: {paper_info['who_umc']}")
            
            paper_details.append(paper_info)
            
        except Exception as e:
            print(f"Error processing paper: {str(e)}")
    
    return paper_details

def plot_combined_analysis(papers, drug_name, qt_data):
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
                line=dict(
                    color=line['color'],
                    dash=line['dash']
                )
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
                marker=dict(
                    size=10,
                    color='red',
                    symbol='circle'
                )
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

    return fig

def process_paper(paper):
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

def save_to_csv(papers, filepath):
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

def extract_qt_hr_values(text):
    """Extract QT and HR values from text."""
    values = {
        'qt_values': [],
        'hr_values': [],
        'bp_values': []
    }
    
    # QT interval patterns
    qt_patterns = [
        r'QT\s*(?:interval)?\s*(?:of|was|is|=|:)?\s*(\d{2,3})(?:\s*(?:ms|msec|milliseconds))?',
        r'QT\s*(?:interval)?\s*(?:prolonged|increased|decreased|measured)\s*(?:to|at)?\s*(\d{2,3})(?:\s*(?:ms|msec|milliseconds))?'
    ]
    
    # Heart rate patterns
    hr_patterns = [
        r'(?:heart rate|HR|pulse|rate)\s*(?:of|was|is|=|:)?\s*(\d{1,3})(?:\s*(?:bpm|beats per minute|beats/min))?',
        r'(?:heart rate|HR|pulse|rate)\s*(?:increased|decreased|measured)\s*(?:to|at)?\s*(\d{1,3})(?:\s*(?:bpm|beats per minute|beats/min))?'
    ]
    
    # Blood pressure patterns
    bp_patterns = [
        r'(?:blood pressure|BP)\s*(?:of|was|is|=|:)?\s*(\d{2,3})/(\d{2,3})(?:\s*(?:mmHg))?',
        r'(?:blood pressure|BP)\s*(?:increased|decreased|measured)\s*(?:to|at)?\s*(\d{2,3})/(\d{2,3})(?:\s*(?:mmHg))?',
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
    
    for pattern in qt_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            try:
                qt = int(match.group(1))
                if 200 <= qt <= 1000:  # Reasonable QT range
                    values['qt_values'].append(qt)
            except ValueError:
                continue

    for pattern in hr_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            try:
                hr = int(match.group(1))
                if 20 <= hr <= 300:  # Reasonable HR range
                    values['hr_values'].append(hr)
            except ValueError:
                continue
    
    # Blood pressure extraction
    for pattern in bp_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            try:
                if len(match.groups()) == 2:
                    sbp = int(match.group(1))
                    dbp = int(match.group(2))
                    if 60 <= sbp <= 250 and 30 <= dbp <= 150:  # Reasonable BP ranges
                        values['bp_values'].append((sbp, dbp))
                else:
                    sbp = int(match.group(1))
                    if 60 <= sbp <= 250:  # Reasonable systolic BP range
                        values['bp_values'].append(sbp)
            except ValueError:
                continue
    
    # Remove duplicates while preserving order
    values['qt_values'] = list(dict.fromkeys(values['qt_values']))
    values['hr_values'] = list(dict.fromkeys(values['hr_values']))
    values['bp_values'] = list(dict.fromkeys(values['bp_values']))
    
    return values

def extract_dose(text):
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

def extract_medical_history(text):
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

def extract_medications(text):
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

def extract_treatment_course(text):
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

def main():
    """Main function to run the TdP risk analysis."""
    import argparse
    
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Analyze TdP risk for a given drug')
    parser.add_argument('drug_name', help='Name of the drug to analyze')
    parser.add_argument('--csv', help='Path to CSV file with QT/HR data', default='qt_data.csv')
    parser.add_argument('--output', help='Output directory', default='output')
    args = parser.parse_args()
    
    drug_name = args.drug_name
    csv_file = args.csv
    output_dir = args.output
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save to docs directory for GitHub Pages
    output_dir = os.path.join(os.path.dirname(__file__), 'docs')
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for hERG and bioavailability data
    print(f"\nSearching for {drug_name} properties...")
    herg_query = f"({drug_name}[Title]) AND (hERG OR IC50 OR potassium channel)"
    bio_query = f"({drug_name}[Title]) AND (bioavailability OR absorption)"
    
    # Initialize properties
    herg_ic50 = None
    bioavailability = None
    
    try:
        # Search PubMed for hERG data
        handle = Entrez.esearch(db="pubmed", term=herg_query)
        record = Entrez.read(handle)
        if record["Count"] != "0":
            for pmid in record["IdList"][:5]:  # Check first 5 papers
                abstract = get_paper_abstract(pmid)
                if abstract:
                    # Look for IC50 values
                    ic50_match = re.search(r'IC50\s*[=:]\s*(\d+\.?\d*)\s*(?:n[Mm]|µ[Mm]|microM)', abstract)
                    if ic50_match:
                        herg_ic50 = float(ic50_match.group(1))
                        break
        
        # Search for bioavailability
        handle = Entrez.esearch(db="pubmed", term=bio_query)
        record = Entrez.read(handle)
        if record["Count"] != "0":
            for pmid in record["IdList"][:5]:
                abstract = get_paper_abstract(pmid)
                if abstract:
                    # Look for bioavailability percentage
                    bio_match = re.search(r'(?:oral\s+)?bioavailability\s*(?:of|was|is)?\s*(\d+\.?\d*)%?', abstract, re.IGNORECASE)
                    if bio_match:
                        bioavailability = float(bio_match.group(1))
                        break
    
    except Exception as e:
        print(f"Error searching for drug properties: {e}")
    
    if herg_ic50:
        print(f"Found hERG IC50: {herg_ic50} nM")
    if bioavailability:
        print(f"Found oral bioavailability: {bioavailability}%")
    
    # Process CSV data
    qt_data = []
    if os.path.exists(csv_file):
        print(f"\nProcessing data from {csv_file}...")
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    hr = float(row['HR'])
                    qt = float(row['QT'])
                    qt_data.append({'hr': hr, 'qt': qt})
                    print(f"Successfully added point: HR={hr}, QT={qt}")
                except (ValueError, KeyError) as e:
                    print(f"Error processing row: {row}, Error: {str(e)}")
        print(f"\nTotal points processed: {len(qt_data)}")
    
    # Search PubMed for TdP cases
    mesh_terms = f"((hERG) OR (QT) OR (QTc) OR (torsad*)) AND ({drug_name})"
    title_abs_terms = f"({drug_name}[Title/Abstract]) AND ((\"torsades de pointes\"[Title/Abstract]) OR (tdp[Title/Abstract]) OR (torsades[Title/Abstract]))"
    
    print(f"\nSearching PubMed using:")
    print(f"MESH terms: {mesh_terms}")
    print(f"Title/Abstract terms: {title_abs_terms}")
    
    papers = search_pubmed(mesh_terms, title_abs_terms)
    print(f"\nFound {len(papers)} matching papers\n")
    
    # Process each paper to extract detailed information
    processed_papers = []
    for paper in papers:
        processed_paper = process_paper(paper)
        
        # Extract QT/HR values from the processed paper
        if processed_paper['QT_values'] != 'N/A' and processed_paper['HR_values'] != 'N/A':
            # Convert string values to lists if they're not already
            qt_values = processed_paper['QT_values'] if isinstance(processed_paper['QT_values'], list) else [processed_paper['QT_values']]
            hr_values = processed_paper['HR_values'] if isinstance(processed_paper['HR_values'], list) else [processed_paper['HR_values']]
            
            # Add each QT/HR pair to qt_data
            for qt, hr in zip(qt_values, hr_values):
                try:
                    qt_val = float(qt)
                    hr_val = float(hr)
                    if 0 < hr_val < 300 and 200 < qt_val < 800:  # Basic validation
                        qt_data.append({'hr': hr_val, 'qt': qt_val})
                        print(f"Added QT/HR point from paper: HR={hr_val}, QT={qt_val}")
                except (ValueError, TypeError):
                    continue
    
    # Save papers to CSV
    csv_path = os.path.join(output_dir, f"{drug_name}_analysis.csv")
    save_to_csv(processed_papers, csv_path)
    print(f"Saved analysis to: {csv_path}")
    
    # Create and save the combined plot
    fig = plot_combined_analysis(processed_papers, drug_name, qt_data)
    html_path = os.path.join(output_dir, f"{drug_name}_analysis.html")
    fig.write_html(html_path)
    print(f"Saved combined analysis to: {html_path}")
    
    # Create index.html that redirects to the analysis
    index_path = os.path.join(output_dir, 'index.html')
    with open(index_path, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <title>QT Nomogram Analysis</title>
    <meta http-equiv="refresh" content="0; url=./{drug_name.lower()}_analysis.html">
</head>
<body>
    <p>Redirecting to <a href="./{drug_name.lower()}_analysis.html">QT Nomogram Analysis</a>...</p>
</body>
</html>
""")

if __name__ == "__main__":
    main()
