import os
import matplotlib.pyplot as plt
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

# Set your email for PubMed API
Entrez.email = "ghhercock@gmail.com"  # Replace with your email

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
    
    return nomogram_line_solid, nomogram_line_dashed

def is_point_above_line(hr, qt):
    # Points with HR > 150 are automatically considered "above" the line
    if hr > 150:
        return True
        
    # Get reference lines
    line_solid, line_dashed = get_reference_lines()
    
    # Convert to numpy arrays for interpolation
    hr_all = [point['hr'] for point in line_solid + line_dashed]
    qt_all = [point['qt'] for point in line_solid + line_dashed]
    
    # Interpolate the line at the given heart rate
    qt_line = np.interp(hr, hr_all, qt_all)
    
    # Point is above or on the line if its QT value is greater than or equal to the line's QT value
    return qt >= qt_line

def plot_qt_nomogram(data_points):
    plt.figure(figsize=(10, 8))  # Approximate 800x550 pixels
    
    # Set margins
    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.15)
    
    # Get reference lines
    line_solid, line_dashed = get_reference_lines()
    
    # Plot grid
    plt.grid(True, linestyle='--', alpha=0.7, dashes=(3, 3))
    
    # Plot reference lines
    hr_solid = [point['hr'] for point in line_solid]
    qt_solid = [point['qt'] for point in line_solid]
    hr_dashed = [point['hr'] for point in line_dashed]
    qt_dashed = [point['qt'] for point in line_dashed]
    
    plt.plot(hr_solid, qt_solid, '-', color='black', linewidth=3)
    plt.plot(hr_dashed, qt_dashed, '--', color='black', linewidth=3, dashes=(5, 5))
    
    # Separate points above and below the line
    above_points = []
    below_points = []
    
    for point in data_points:
        if is_point_above_line(point['hr'], point['qt']):
            above_points.append(point)
        else:
            below_points.append(point)
    
    # Plot points above the line in red
    if above_points:
        hr_above = [point['hr'] for point in above_points]
        qt_above = [point['qt'] for point in above_points]
        plt.scatter(hr_above, qt_above, color='red', s=30)
    
    # Plot points below the line in black
    if below_points:
        hr_below = [point['hr'] for point in below_points]
        qt_below = [point['qt'] for point in below_points]
        plt.scatter(hr_below, qt_below, color='black', s=30)
    
    # Set axis limits and ticks
    plt.xlim(0, 260)  # Increased to show HR values up to 250
    plt.ylim(200, 900)
    
    plt.xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260])
    plt.yticks([200, 300, 400, 500, 600, 700, 800, 900])
    
    # Add labels with bold font
    plt.xlabel('Heart Rate (bpm)', fontweight='bold', labelpad=15)
    plt.ylabel('QT Interval (ms)', fontweight='bold')
    
    # Style the axes
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    
    # Save the plot with high DPI for clarity
    plt.savefig('qt_nomogram.png', dpi=300, bbox_inches='tight')
    plt.close()

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
    dose_pattern = r'(\d+(?:\.\d+)?)\s*(?:-\s*\d+(?:\.\d+)?)?\s*mg'
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
    fig, ax = plt.subplots(figsize=(20, len(data) + 1))
    
    # Hide axes
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    table = ax.table(cellText=data,
                    colLabels=columns,
                    cellLoc='left',
                    loc='center',
                    colWidths=[0.3] + [0.07] * (len(columns)-1))
    
    # Adjust font size and row height
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    
    # Add title
    plt.title(f'{drug_name} Case Reports Analysis', pad=20)
    
    # Save figure
    plt.savefig(f'{drug_name.lower()}_analysis_table.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{drug_name.lower()}_analysis_table.pdf', bbox_inches='tight')
    plt.close()

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
                            for aid in record['AID']:
                                if '[doi]' in aid:
                                    doi = aid.replace(' [doi]', '')
                                    break
                        
                        paper = {
                            'title': record.get('TI', ''),
                            'abstract': record.get('AB', ''),
                            'authors': ', '.join(record.get('AU', [])),
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
            qt_hr_values = extract_qt_hr_values(text)
            
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
                'QT (ms)': qt_hr_values['uncorrected_qt'] or 'N/A',
                'QTc (ms)': qt_hr_values['qtc'] or 'N/A',
                'QTB (ms)': qt_hr_values['qtb'] or 'N/A',
                'QTF (ms)': qt_hr_values['qtf'] or 'N/A',
                'HR (bpm)': qt_hr_values['hr'] or 'N/A',
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
            paper_info['naranjo_score'] = calculate_naranjo_score(paper_info)
            paper_info['tisdale_score'] = calculate_tisdale_score(paper_info)
            paper_info['who_umc'] = calculate_who_umc_score(paper_info)
            
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
    # Create figure with two subplots - adjust height based on number of papers
    fig = plt.figure(figsize=(20, 10 + len(papers)))
    
    # Add QT nomogram subplot
    ax1 = plt.subplot2grid((2, 1), (0, 0))
    
    # Get reference lines
    line_solid, line_dashed = get_reference_lines()
    
    # Plot grid
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Plot reference lines
    hr_solid = [point['hr'] for point in line_solid]
    qt_solid = [point['qt'] for point in line_solid]
    hr_dashed = [point['hr'] for point in line_dashed]
    qt_dashed = [point['qt'] for point in line_dashed]
    
    ax1.plot(hr_solid, qt_solid, '-', color='black', linewidth=2, label='Reference Line')
    ax1.plot(hr_dashed, qt_dashed, '--', color='black', linewidth=2)
    
    # Plot data points
    if qt_data:
        hr_values = [point['hr'] for point in qt_data]
        qt_values = [point['qt'] for point in qt_data]
        ax1.scatter(hr_values, qt_values, color='red', s=50, label='Data Points')
    
    ax1.set_xlabel('Heart Rate (bpm)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('QT Interval (ms)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{drug_name} QT Nomogram Analysis', fontsize=14, fontweight='bold', pad=20)
    ax1.legend(fontsize=10)
    
    # Style the axes
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(labelsize=10)
    
    # Set axis limits with padding
    ax1.set_xlim(0, 160)
    ax1.set_ylim(200, 800)
    
    # Extract table data
    data = []
    columns = ['Title', 'Year', 'DOI', 'Authors', 'Age', 'Sex', 'Dose (mg)', 
              'QT (ms)', 'QTc (ms)', 'QTB (ms)', 'QTF (ms)', 'HR (bpm)', 
              'TdP?', 'BP (mmHg)', 'Medical History', 'Medication History',
              'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC']
    
    for paper in papers:
        text = paper['title'] + ' ' + paper['abstract']
        
        # Extract all information
        dose = extract_dose(text)
        qt_hr_values = extract_qt_hr_values(text)
        
        # Extract age and sex with capture groups
        age_match = re.search(r'(\d+)[\s-]*(year|yr)', text, re.IGNORECASE)
        sex_match = re.search(r'(male|female|man|woman|boy|girl)', text, re.IGNORECASE)
        
        row = [
            paper['title'],
            paper['year'],
            paper.get('doi', 'N/A'),
            paper['authors'],
            age_match.group(1) if age_match else 'N/A',
            sex_match.group(1).title() if sex_match else 'N/A',
            str(dose) if dose else 'N/A',
            qt_hr_values['uncorrected_qt'] or 'N/A',
            qt_hr_values['qtc'] or 'N/A',
            qt_hr_values['qtb'] or 'N/A',
            qt_hr_values['qtf'] or 'N/A',
            qt_hr_values['hr'] or 'N/A',
            'Yes' if 'torsade' in text.lower() else 'N/A',
            qt_hr_values['bp'] or 'N/A',
            extract_medical_history(text),
            extract_medications(text),
            extract_treatment_course(text),
            extract_outcome(text),
            str(paper['naranjo_score']),
            str(paper['tisdale_score']),
            paper['who_umc']
        ]
        data.append(row)
    
    # Add table subplot with auto-sized columns
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    ax2.axis('tight')
    ax2.axis('off')
    
    # Create table with auto column widths
    table = ax2.table(cellText=data,
                     colLabels=columns,
                     cellLoc='left',
                     loc='center')
    
    # Auto-size columns based on content
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    
    # Calculate column widths based on content
    widths = []
    for col in range(len(columns)):
        column_contents = [str(row[col]) for row in data]
        column_contents.append(columns[col])
        max_width = max(len(str(content)) for content in column_contents)
        widths.append(max_width)
    
    # Normalize widths to sum to ~1
    total_width = sum(widths)
    col_widths = [w/total_width for w in widths]
    
    # Apply column widths and style cells
    for (row, col), cell in table.get_celld().items():
        cell.set_width(col_widths[col])
        if row == 0:
            cell.set_text_props(rotation=45)
            cell.set_height(0.15)
            cell.set_text_props(weight='bold')
        else:
            cell.set_height(0.1)
            # Add word wrapping for text columns
            if col in [0, 3, 14, 15, 16]:  # Title, Authors, and text fields
                cell._text.set_wrap(True)
    
    # Scale table to fit
    table.scale(1, 2)
    
    plt.title(f'{drug_name} Case Reports Analysis', pad=20)
    plt.tight_layout()
    
    # Save with high DPI
    plt.savefig(f'{drug_name.lower()}_full_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{drug_name.lower()}_full_analysis.pdf', bbox_inches='tight')
    plt.close()

def extract_medical_history(text):
    """Extract medical history from text."""
    conditions = []
    
    # Common conditions and their variations
    condition_patterns = [
        # Cardiovascular conditions
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
        r'((?:prolonged|long)\s+QT(?:\s+syndrome)?)',
        # Respiratory conditions
        r'(asthma|copd|bronchitis)',
        r'(respiratory\s+(?:failure|distress|infection))',
        r'(pneumonia)',
        # Neurological conditions
        r'(seizure|epilepsy)',
        r'(stroke|tia)',
        r'(parkinson)',
        # Psychiatric conditions
        r'(depression)',
        r'(anxiety)',
        r'(bipolar\s+disorder)',
        r'(schizophrenia)',
        # Other conditions
        r'(diabetes(?:\s+mellitus)?)',
        r'(renal\s+(?:failure|insufficiency))',
        r'(liver\s+(?:failure|disease|cirrhosis))',
        r'(cancer|malignancy)',
        r'(sepsis|infection)',
        r'(covid-19|coronavirus)',
        # Lab values
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
        # Administration patterns
        r'(?:taking|on|receiving|treated with|therapy with)\s+([\w\s-]+(?:blocker|inhibitor|antagonist|agonist|sartan|pril|statin))',
        r'(?:prescribed|administered|given)\s+([\w\s-]+(?:mg|mcg|g|ml|units))',
        r'(\w+)\s+(?:\d+\s*(?:mg|mcg|g|ml|units))',
        # Drug classes
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
        # Common cardiovascular drugs
        r'\b(warfarin|heparin|aspirin|clopidogrel)',
        r'\b(metoprolol|atenolol|carvedilol|bisoprolol)',
        r'\b(lisinopril|enalapril|ramipril)',
        r'\b(amlodipine|diltiazem|verapamil)',
        r'\b(furosemide|bumetanide|spironolactone)',
        r'\b(amiodarone|sotalol|flecainide)',
        r'\b(digoxin)',
        # Common antibiotics
        r'\b(azithromycin|clarithromycin|erythromycin)',
        r'\b(ciprofloxacin|levofloxacin|moxifloxacin)',
        r'\b(penicillin|amoxicillin|ampicillin)',
        r'\b(ceftriaxone|cefuroxime|cephalexin)',
        # Common psychiatric drugs
        r'\b(fluoxetine|sertraline|paroxetine|citalopram|escitalopram)',
        r'\b(quetiapine|olanzapine|risperidone|haloperidol)',
        r'\b(lithium|valproate|carbamazepine)',
        # Other common drugs
        r'\b(metformin|insulin)',
        r'\b(hydroxychloroquine|chloroquine)',
        r'\b(ivermectin|remdesivir)',
        # Drug combinations
        r'(\w+(?:\s*[+]\s*\w+)+)',
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
        # General treatment patterns
        r'(?:treated|managed)\s+with\s+([\w\s-]+)',
        r'(?:administered|given|received)\s+([\w\s-]+)',
        r'(?:therapy|treatment)\s+with\s+([\w\s-]+)',
        # Specific interventions
        r'(defibrillation|cardioversion|pacing)',
        r'(magnesium\s+(?:infusion|supplementation))',
        r'(potassium\s+(?:infusion|supplementation))',
        r'(activated\s+charcoal)',
        r'(intubation|ventilation)',
        r'(cpr|cardiopulmonary\s+resuscitation)',
        r'(hemodialysis|dialysis)',
        r'(ecg\s+monitoring)',
        r'(telemetry)',
        # Drug adjustments
        r'(discontinuation|withdrawal|stopped)\s+of\s+([\w\s-]+)',
        r'(dose\s+reduction|decreased|increased)\s+of\s+([\w\s-]+)',
        r'(monitoring|observation)',
        # Specific treatments
        r'(iv\s+fluids?)',
        r'(oxygen\s+therapy)',
        r'(mechanical\s+ventilation)',
        r'(vasopressors?)',
        r'(inotropes?)',
        r'(antibiotics?)',
        r'(antiarrhythmics?)',
        # Procedures
        r'(surgery|operation)',
        r'(catheterization)',
        r'(angiography)',
        r'(stent(?:ing)?)',
        r'(bypass)',
        # Supportive care
        r'(supportive\s+care)',
        r'(symptomatic\s+treatment)',
        r'(palliative\s+care)',
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

def extract_outcome(text):
    """Extract patient outcome from text."""
    # Outcome patterns
    outcome_patterns = [
        # Death patterns
        (r'(?:patient|subject)\s+(?:died|expired|deceased|death)', 'Death'),
        (r'(?:fatal|lethal)\s+outcome', 'Death'),
        # Recovery patterns
        (r'(?:complete|full)\s+recovery', 'Complete recovery'),
        (r'discharged\s+(?:home|from\s+hospital)', 'Discharged'),
        (r'returned\s+to\s+baseline', 'Returned to baseline'),
        # Partial recovery
        (r'partial\s+recovery', 'Partial recovery'),
        (r'improved\s+(?:with|after)', 'Improved'),
        # Ongoing/chronic
        (r'chronic|ongoing|persistent', 'Ongoing symptoms'),
        # Transfer
        (r'transferred\s+to', 'Transferred'),
        # Default pattern
        (r'(?:survived|alive)', 'Survived'),
    ]
    
    for pattern, outcome in outcome_patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return outcome
    
    return 'Not reported'

def extract_qt_hr_values(text):
    """Extract QT and HR values using comprehensive patterns."""
    values = {
        'uncorrected_qt': None,
        'qtc': None,
        'qtb': None,
        'qtf': None,
        'hr': None,
        'bp': None
    }
    
    # QT patterns
    qt_patterns = [
        # General QT patterns
        r'(?:QT|qt)\s*(?:interval)?\s*(?:of|was|is|=|:)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
        r'(?:QT|qt)\s*(?:interval)?\s*(?:prolonged|increased|extended)\s*to\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
        # QTc patterns
        r'(?:QTc|qtc|corrected\s*QT)\s*(?:interval)?\s*(?:of|was|is|=|:)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
        r'(?:QTc|qtc|corrected\s*QT)\s*(?:interval)?\s*(?:prolonged|increased|extended)\s*to\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
        # QTB patterns
        r'(?:QTb|qtb|Bazett)\s*(?:interval|correction)?\s*(?:of|was|is|=|:)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
        # QTF patterns
        r'(?:QTf|qtf|Fridericia)\s*(?:interval|correction)?\s*(?:of|was|is|=|:)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:ms|msec|milliseconds)',
    ]
    
    # HR patterns
    hr_patterns = [
        # Direct HR mentions
        r'(?:HR|heart\s*rate|ventricular\s*rate|pulse(?:\s*rate)?)\s*(?:of|was|is|=|:)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:bpm|beats\s*(?:per|/)\s*min|b\.p\.m\.|min-1)',
        r'(?:HR|heart\s*rate|ventricular\s*rate|pulse(?:\s*rate)?)\s*(?:increased|decreased)\s*to\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:bpm|beats\s*(?:per|/)\s*min|b\.p\.m\.|min-1)',
        # Bradycardia/Tachycardia mentions
        r'(?:bradycardia|tachycardia)\s*(?:at|of|with|showing)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:bpm|beats\s*(?:per|/)\s*min|b\.p\.m\.|min-1)',
        # General rate mentions
        r'rate\s*(?:of|at|showing)\s*(\d+)(?:\.|,)?(?:\d+)?\s*(?:bpm|beats\s*(?:per|/)\s*min|b\.p\.m\.|min-1)',
        # Numeric only with units
        r'(\d+)(?:\.|,)?(?:\d+)?\s*(?:bpm|beats\s*(?:per|/)\s*min|b\.p\.m\.|min-1)',
    ]
    
    # BP patterns
    bp_patterns = [
        r'(?:BP|blood\s*pressure)\s*(?:of|was|is|=|:)\s*(\d+)/(\d+)',
        r'(?:systolic|SBP)\s*(?:of|was|is|=|:)\s*(\d+)',
        r'(?:diastolic|DBP)\s*(?:of|was|is|=|:)\s*(\d+)'
    ]
    
    # Extract QT values
    for pattern in qt_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            if 'QTc' in pattern or 'qtc' in pattern or 'corrected' in pattern:
                values['qtc'] = int(match.group(1))
            elif 'QTb' in pattern or 'qtb' in pattern or 'Bazett' in pattern:
                values['qtb'] = int(match.group(1))
            elif 'QTf' in pattern or 'qtf' in pattern or 'Fridericia' in pattern:
                values['qtf'] = int(match.group(1))
            else:
                values['uncorrected_qt'] = int(match.group(1))
    
    # Extract HR value
    for pattern in hr_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            try:
                hr = float(match.group(1))
                if 20 <= hr <= 300:  # Reasonable HR range
                    values['hr'] = hr
                    break
            except ValueError:
                continue
    
    # Extract BP values
    for pattern in bp_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            if len(match.groups()) == 2:
                values['bp'] = f"{match.group(1)}/{match.group(2)}"
            else:
                values['bp'] = match.group(1)
            break
    
    return values

def calculate_naranjo_score(text):
    """Calculate Naranjo score for adverse drug reaction causality."""
    score = 0
    
    # Previous reports
    if re.search(r'previous.*report|case.*report|literature|published', text, re.I):
        score += 1
    
    # Temporal relationship
    if re.search(r'after|following|subsequent|began|started', text, re.I):
        score += 2
    
    # Improvement after discontinuation
    if re.search(r'discontinu|withdraw|stop|ceased', text, re.I):
        if re.search(r'improv|resolv|recover|normal', text, re.I):
            score += 1
    
    # Recurrence on rechallenge
    if re.search(r'rechallenge|reexposure', text, re.I):
        score += 2
    
    # Alternative causes
    if not re.search(r'alternative|other cause', text, re.I):
        score += 2
    
    # Dose-response relationship
    if re.search(r'dose.*relationship|concentration|level|plasma', text, re.I):
        score += 1
    
    # Objective evidence
    if re.search(r'ECG|EKG|monitor|measure|recorded', text, re.I):
        score += 1
    
    return score

def calculate_tisdale_score(text):
    """Calculate Tisdale score for QT prolongation risk."""
    score = 0
    
    # Age ≥68 years
    age_match = re.search(r'(\d+)[\s-]*(?:year|yr)', text, re.I)
    if age_match and int(age_match.group(1)) >= 68:
        score += 1
    
    # Female sex
    if re.search(r'female|woman', text, re.I):
        score += 1
    
    # Loop diuretic
    if re.search(r'furosemide|bumetanide|torsemide', text, re.I):
        score += 1
    
    # Serum K+ ≤3.5
    k_match = re.search(r'potassium.*?(\d+\.?\d*)', text, re.I)
    if k_match and float(k_match.group(1)) <= 3.5:
        score += 2
    
    # QTc on admission ≥450 ms
    qtc_match = re.search(r'QTc.*?(\d+)', text)
    if qtc_match and int(qtc_match.group(1)) >= 450:
        score += 2
    
    # Acute MI
    if re.search(r'myocardial infarction|MI|heart attack', text, re.I):
        score += 2
    
    # Sepsis
    if re.search(r'sepsis|septic', text, re.I):
        score += 3
    
    # Heart failure
    if re.search(r'heart failure|CHF|cardiac failure', text, re.I):
        score += 3
    
    # One QTc-prolonging drug
    if re.search(r'antiarrhythmic|antipsychotic|antibiotic|antidepressant', text, re.I):
        score += 3
    
    # ≥2 QTc-prolonging drugs
    drug_count = len(re.findall(r'antiarrhythmic|antipsychotic|antibiotic|antidepressant', text, re.I))
    if drug_count >= 2:
        score += 3
    
    return score

def assess_who_umc(text):
    """Assess WHO-UMC causality category."""
    # Certain
    if (re.search(r'rechallenge|reexposure', text, re.I) and 
        re.search(r'positive|confirm|definite', text, re.I)):
        return "Certain"
    
    # Probable/Likely
    if (re.search(r'temporal|relationship|association', text, re.I) and 
        re.search(r'improve|resolve|recover', text, re.I) and 
        not re.search(r'alternative|other cause', text, re.I)):
        return "Probable/Likely"
    
    # Possible
    if (re.search(r'temporal|relationship|association', text, re.I) and 
        re.search(r'could|may|possible', text, re.I)):
        return "Possible"
    
    # Unlikely
    if re.search(r'unlikely|doubtful|improbable', text, re.I):
        return "Unlikely"
    
    # Default to Possible if no clear category
    return "Possible"

def download_pdf_from_doi(doi, output_dir):
    """Try to download PDF from Sci-Hub or other open access sources."""
    try:
        # First try unpaywall API
        unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email=your.email@example.com"
        response = requests.get(unpaywall_url)
        if response.status_code == 200:
            data = response.json()
            if data.get('is_oa') and data.get('best_oa_location', {}).get('url_for_pdf'):
                pdf_url = data['best_oa_location']['url_for_pdf']
                response = requests.get(pdf_url)
                if response.status_code == 200:
                    filename = f"{output_dir}/{doi.replace('/', '_')}.pdf"
                    with open(filename, 'wb') as f:
                        f.write(response.content)
                    print(f"Successfully downloaded: {filename}")
                    return filename
        
        print(f"Could not find free PDF for DOI: {doi}")
        return None
    
    except Exception as e:
        print(f"Error downloading PDF: {str(e)}")
        return None

def search_and_download_papers(drug_name, output_dir):
    """Search for papers and try to download their PDFs."""
    papers = search_pubmed_papers(drug_name)
    
    if not papers:
        return []
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save paper details to CSV
    csv_path = os.path.join(output_dir, f"{drug_name}_papers.csv")
    save_to_csv(papers, csv_path)
    print(f"\nSaved paper details to: {csv_path}")
    
    # Try to download PDFs
    downloaded_pdfs = []
    for paper in papers:
        if 'DOI' in paper:
            pdf_path = download_pdf_from_doi(paper['DOI'], output_dir)
            if pdf_path:
                downloaded_pdfs.append(pdf_path)
        time.sleep(1)  # Be nice to the servers
    
    return downloaded_pdfs

def plot_combined_analysis(papers, drug_name, qt_data):
    """Create a combined figure with QT nomogram and analysis table."""
    # Create figure with two subplots - adjust height based on number of papers
    fig = plt.figure(figsize=(20, 10 + len(papers)))
    
    # Add QT nomogram subplot
    ax1 = plt.subplot2grid((2, 1), (0, 0))
    
    # Get reference lines
    line_solid, line_dashed = get_reference_lines()
    
    # Plot grid
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Plot reference lines
    hr_solid = [point['hr'] for point in line_solid]
    qt_solid = [point['qt'] for point in line_solid]
    hr_dashed = [point['hr'] for point in line_dashed]
    qt_dashed = [point['qt'] for point in line_dashed]
    
    ax1.plot(hr_solid, qt_solid, '-', color='black', linewidth=2, label='Reference Line')
    ax1.plot(hr_dashed, qt_dashed, '--', color='black', linewidth=2)
    
    # Plot data points
    if qt_data:
        hr_values = [point['hr'] for point in qt_data]
        qt_values = [point['qt'] for point in qt_data]
        ax1.scatter(hr_values, qt_values, color='red', s=50, label='Data Points')
    
    ax1.set_xlabel('Heart Rate (bpm)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('QT Interval (ms)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{drug_name} QT Nomogram Analysis', fontsize=14, fontweight='bold', pad=20)
    ax1.legend(fontsize=10)
    
    # Style the axes
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(labelsize=10)
    
    # Set axis limits with padding
    ax1.set_xlim(0, 160)
    ax1.set_ylim(200, 800)
    
    # Extract table data
    data = []
    columns = ['Title', 'Year', 'DOI', 'Authors', 'Age', 'Sex', 'Dose (mg)', 
              'QT (ms)', 'QTc (ms)', 'QTB (ms)', 'QTF (ms)', 'HR (bpm)', 
              'TdP?', 'BP (mmHg)', 'Medical History', 'Medication History',
              'Course of Treatment', 'Outcome', 'Naranjo', 'Tisdale', 'WHO-UMC']
    
    for paper in papers:
        text = paper['title'] + ' ' + paper['abstract']
        
        # Extract all information
        dose = extract_dose(text)
        qt_hr_values = extract_qt_hr_values(text)
        
        # Extract age and sex with capture groups
        age_match = re.search(r'(\d+)[\s-]*(year|yr)', text, re.IGNORECASE)
        sex_match = re.search(r'(male|female|man|woman|boy|girl)', text, re.IGNORECASE)
        
        row = [
            paper['title'],
            paper['year'],
            paper.get('doi', 'N/A'),
            paper['authors'],
            age_match.group(1) if age_match else 'N/A',
            sex_match.group(1).title() if sex_match else 'N/A',
            str(dose) if dose else 'N/A',
            qt_hr_values['uncorrected_qt'] or 'N/A',
            qt_hr_values['qtc'] or 'N/A',
            qt_hr_values['qtb'] or 'N/A',
            qt_hr_values['qtf'] or 'N/A',
            qt_hr_values['hr'] or 'N/A',
            'Yes' if 'torsade' in text.lower() else 'N/A',
            qt_hr_values['bp'] or 'N/A',
            extract_medical_history(text),
            extract_medications(text),
            extract_treatment_course(text),
            extract_outcome(text),
            str(paper['naranjo_score']),
            str(paper['tisdale_score']),
            paper['who_umc']
        ]
        data.append(row)
    
    # Add table subplot with auto-sized columns
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    ax2.axis('tight')
    ax2.axis('off')
    
    # Create table with auto column widths
    table = ax2.table(cellText=data,
                     colLabels=columns,
                     cellLoc='left',
                     loc='center')
    
    # Auto-size columns based on content
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    
    # Calculate column widths based on content
    widths = []
    for col in range(len(columns)):
        column_contents = [str(row[col]) for row in data]
        column_contents.append(columns[col])
        max_width = max(len(str(content)) for content in column_contents)
        widths.append(max_width)
    
    # Normalize widths to sum to ~1
    total_width = sum(widths)
    col_widths = [w/total_width for w in widths]
    
    # Apply column widths and style cells
    for (row, col), cell in table.get_celld().items():
        cell.set_width(col_widths[col])
        if row == 0:
            cell.set_text_props(rotation=45)
            cell.set_height(0.15)
            cell.set_text_props(weight='bold')
        else:
            cell.set_height(0.1)
            # Add word wrapping for text columns
            if col in [0, 3, 14, 15, 16]:  # Title, Authors, and text fields
                cell._text.set_wrap(True)
    
    # Scale table to fit
    table.scale(1, 2)
    
    plt.title(f'{drug_name} Case Reports Analysis', pad=20)
    plt.tight_layout()
    
    # Save with high DPI
    plt.savefig(f'{drug_name.lower()}_full_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{drug_name.lower()}_full_analysis.pdf', bbox_inches='tight')
    plt.close()

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
                    print(f"Error processing row: {e}")
        print(f"\nTotal points processed: {len(qt_data)}")
    
    # Search PubMed for TdP cases
    mesh_terms = f"((hERG) OR (QT) OR (QTc) OR (torsad*)) AND ({drug_name})"
    title_abs_terms = f"({drug_name}[Title/Abstract]) AND ((\"torsades de pointes\"[Title/Abstract]) OR (tdp[Title/Abstract]) OR (torsades[Title/Abstract]))"
    
    print(f"\nSearching PubMed using:")
    print(f"MESH terms: {mesh_terms}")
    print(f"Title/Abstract terms: {title_abs_terms}")
    
    papers = search_pubmed(mesh_terms, title_abs_terms)
    print(f"\nFound {len(papers)} matching papers\n")
    
    # Extract QT/HR data from papers
    for paper in papers:
        text = paper['title'] + ' ' + paper['abstract']
        qt_hr_values = extract_qt_hr_values(text)
        if qt_hr_values['hr'] and qt_hr_values['uncorrected_qt']:
            qt_data.append({
                'hr': float(qt_hr_values['hr']),
                'qt': float(qt_hr_values['uncorrected_qt'])
            })
    
    # Save papers to CSV
    csv_path = os.path.join(output_dir, f"{drug_name}_analysis.csv")
    save_to_csv(papers, csv_path)
    print(f"Saved analysis to: {csv_path}")
    
    # Create and save the nomogram plot
    plot_path = os.path.join(output_dir, f"{drug_name}_nomogram.png")
    create_nomogram_plot(drug_name, qt_data, plot_path, herg_ic50, bioavailability)
    print(f"Saved nomogram plot to: {plot_path}")

def create_nomogram_plot(drug_name, qt_data, output_path, herg_ic50=None, bioavailability=None):
    """Create QT nomogram plot with data points."""
    plt.figure(figsize=(12, 8))
    
    # Get reference lines
    line_solid, line_dashed = get_reference_lines()
    
    # Plot grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Plot reference lines
    hr_solid = [point['hr'] for point in line_solid]
    qt_solid = [point['qt'] for point in line_solid]
    hr_dashed = [point['hr'] for point in line_dashed]
    qt_dashed = [point['qt'] for point in line_dashed]
    
    plt.plot(hr_solid, qt_solid, '-', color='black', linewidth=2, label='Reference Line')
    plt.plot(hr_dashed, qt_dashed, '--', color='black', linewidth=2)
    
    # Plot data points if available
    if qt_data:
        hr_values = [point['hr'] for point in qt_data]
        qt_values = [point['qt'] for point in qt_data]
        plt.scatter(hr_values, qt_values, color='red', s=50, label='Data Points')
    
    plt.xlabel('Heart Rate (bpm)', fontsize=12, fontweight='bold')
    plt.ylabel('QT Interval (ms)', fontsize=12, fontweight='bold')
    
    # Add title with drug properties if available
    title = f'{drug_name} QT Nomogram Analysis'
    if herg_ic50 or bioavailability:
        props = []
        if herg_ic50:
            props.append(f'hERG IC50: {herg_ic50} nM')
        if bioavailability:
            props.append(f'Bioavailability: {bioavailability}%')
        title += f'\n({", ".join(props)})'
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.legend(fontsize=10)
    
    # Style the axes
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tick_params(labelsize=10)
    
    # Set axis limits with padding
    plt.xlim(0, 160)
    plt.ylim(200, 800)
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
