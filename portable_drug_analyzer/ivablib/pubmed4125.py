import os
import logging
import pandas as pd
import requests
import urllib3
from bs4 import BeautifulSoup
from Bio import Entrez
from Bio import Medline
from typing import List, Dict, Optional
import sys
import time
import re
from urllib.parse import quote
import PyPDF2
import difflib
import csv

# Suppress SSL warnings
urllib3.disable_warnings()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email for NCBI
Entrez.email = "your@email.com"

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'}

def get_from_scihub(url: str, output_dir: str, filename: str) -> Optional[str]:
    """Try to download a paper from Sci-Hub"""
    try:
        # Clean up DOI URL
        clean_url = url.replace('https://doi.org/', '')
        if not clean_url.startswith('http'):
            clean_url = f"https:/doi.org/{clean_url}"
        
        # Try different Sci-Hub domains
        domains = [
            "https://www.sci-hub.se",
            "https://www.sci-hub.ee",
            "https://2024.sci-hub.se",
            "https://dacemirror.sci-hub.se"
        ]
        
        for domain in domains:
            try:
                sci_hub_url = f"{domain}/{clean_url}"
                logger.info(f"Trying URL: {sci_hub_url}")
                
                # Get the Sci-Hub page
                response = requests.get(sci_hub_url, headers=HEADERS, timeout=30)
                if response.status_code != 200:
                    continue
                    
                # Parse the page
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Try to find PDF URL in embed tag
                embed = soup.find('embed')
                if embed and 'src' in embed.attrs:
                    pdf_url = embed['src']
                    if pdf_url.startswith('//'):
                        pdf_url = 'https:' + pdf_url
                    elif not pdf_url.startswith('http'):
                        # Add domain to relative URLs
                        if pdf_url.startswith('/'):
                            pdf_url = domain + pdf_url
                        else:
                            # Special handling for dacemirror
                            if 'dacemirror' in pdf_url:
                                pdf_url = 'https://' + pdf_url
                            else:
                                pdf_url = domain + '/' + pdf_url
                    
                    # Remove URL parameters
                    pdf_url = pdf_url.split('#')[0]
                    logger.info(f"Found PDF URL (embed): {pdf_url}")
                    
                    # Download the PDF
                    pdf_path = os.path.join(output_dir, filename)
                    pdf_response = requests.get(pdf_url, headers=HEADERS, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        if not pdf_response.content.startswith(b'%PDF'):
                            logger.warning(f"Downloaded file is not a valid PDF from {pdf_url}")
                            continue
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
                        logger.info(f"Saved PDF to {pdf_path}")
                        return pdf_path
                
                # Try to find PDF URL in iframe
                iframe = soup.find('iframe')
                if iframe and 'src' in iframe.attrs:
                    pdf_url = iframe['src']
                    if pdf_url.startswith('//'):
                        pdf_url = 'https:' + pdf_url
                    elif not pdf_url.startswith('http'):
                        # Add domain to relative URLs
                        if pdf_url.startswith('/'):
                            pdf_url = domain + pdf_url
                        else:
                            # Special handling for dacemirror
                            if 'dacemirror' in pdf_url:
                                pdf_url = 'https://' + pdf_url
                            else:
                                pdf_url = domain + '/' + pdf_url
                    
                    # Remove URL parameters
                    pdf_url = pdf_url.split('#')[0]
                    logger.info(f"Found PDF URL (iframe): {pdf_url}")
                    
                    # Download the PDF
                    pdf_path = os.path.join(output_dir, filename)
                    pdf_response = requests.get(pdf_url, headers=HEADERS, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        if not pdf_response.content.startswith(b'%PDF'):
                            logger.warning(f"Downloaded file is not a valid PDF from {pdf_url}")
                            continue
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
                        logger.info(f"Saved PDF to {pdf_path}")
                        return pdf_path
                            
            except Exception as e:
                logger.warning(f"Error processing {domain}: {e}")
                continue
                
        return None
        
    except Exception as e:
        logger.error(f"Error downloading from Sci-Hub: {e}")
        return None

def get_doi_from_crossref(title: str) -> Optional[str]:
    """Try to find DOI using CrossRef API with improved title matching"""
    try:
        logger.info(f"Searching CrossRef for: {title}")
        encoded_title = quote(title)
        url = f"https://api.crossref.org/works?query.title={encoded_title}&rows=5&select=DOI,title,score"
        
        headers = {
            'User-Agent': 'IvabradineResearch/1.0 (mailto:your@email.com)',
            'Accept': 'application/json'
        }
        
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()  # Raise exception for non-200 status codes
        
        data = response.json()
        items = data.get('message', {}).get('items', [])
        
        if items:
            # Calculate similarity scores between search title and results
            title_lower = title.lower()
            best_match = None
            best_score = 0
            
            for item in items:
                item_title = item.get('title', [''])[0].lower()
                # Use both CrossRef's score and string similarity
                crossref_score = item.get('score', 0)
                similarity = difflib.SequenceMatcher(None, title_lower, item_title).ratio()
                combined_score = crossref_score * similarity
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_match = item
            
            # Only return DOI if the match is good enough
            if best_match and best_score > 0.5:
                doi = best_match.get('DOI')
                if doi:
                    logger.info(f"Found DOI from CrossRef: {doi} (score: {best_score:.2f})")
                    return doi
            else:
                logger.warning(f"No good match found in CrossRef results (best score: {best_score:.2f})")
                    
    except requests.exceptions.RequestException as e:
        logger.error(f"CrossRef request error for {title}: {e}")
    except ValueError as e:
        logger.error(f"CrossRef JSON parsing error for {title}: {e}")
    except Exception as e:
        logger.error(f"Unexpected CrossRef error for {title}: {e}")
    
    return None

def get_paper_pdf(paper: Dict, output_dir: str) -> Optional[str]:
    """Get PDF for a paper, trying multiple sources"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate safe filename from title
        title = paper.get('title', '')
        safe_title = re.sub(r'[^\w\s-]', '', title)[:100]  # Limit length and remove special chars
        pdf_path = os.path.join(output_dir, f"{safe_title}.pdf")
        
        # Skip if PDF already exists
        if os.path.exists(pdf_path):
            logger.info(f"PDF already exists: {pdf_path}")
            return pdf_path
        
        # Try PMC first if available
        pmc_id = paper.get('pmc')
        if pmc_id:
            logger.info(f"Found PMC ID: {pmc_id}")
            pdf_url = get_pmc_pdf_url(pmc_id)
            if pdf_url:
                if download_pdf(pdf_url, pdf_path):
                    return pdf_path
        
        # Try DOI next
        doi = paper.get('doi')
        if not doi:
            # Try to get DOI from CrossRef
            doi = get_doi_from_crossref(title)
            if doi:
                paper['doi'] = doi
        
        if doi:
            # Try Sci-Hub
            logger.info(f"Trying Sci-Hub for DOI: {doi}")
            pdf_path = get_from_scihub(doi, output_dir, f"{paper.get('pmid', '')}.pdf")
            if pdf_path:
                return pdf_path
        
        logger.warning(f"Could not find PDF for paper: {title}")
        return None
        
    except Exception as e:
        logger.error(f"Error getting PDF for paper {paper.get('title', '')}: {e}")
        return None

def extract_text_from_pdf(pdf_path: str) -> Optional[str]:
    """Extract text from PDF file"""
    try:
        with open(pdf_path, 'rb') as file:
            # Create PDF reader object
            pdf_reader = PyPDF2.PdfReader(file)
            
            # Extract text from all pages
            text = ""
            for page in pdf_reader.pages:
                text += page.extract_text() + "\n"
                
            return text
    except Exception as e:
        logger.warning(f"Failed to extract text from PDF {pdf_path}: {e}")
        return None

def get_drug_names(drug_name: str) -> List[str]:
    """Get all related drug names from the drugnames.csv file"""
    try:
        # Look for CSV in data directory
        csv_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'drugnames.csv')
        if not os.path.exists(csv_path):
            logger.warning(f"drugnames.csv not found at {csv_path}")
            return [drug_name]
            
        # Read CSV and find matches
        df = pd.read_csv(csv_path)
        matches = df[df['DrugName'].str.lower().str.contains(drug_name.lower(), na=False)]
        
        if matches.empty:
            logger.warning(f"No matches found for {drug_name} in drugnames.csv")
            return [drug_name]
            
        # Get both drug name and active ingredient
        drug_names = set()
        for _, row in matches.iterrows():
            drug_names.add(row['DrugName'])
            if pd.notna(row['ActiveIngredient']):
                drug_names.add(row['ActiveIngredient'])
                
        logger.info(f"Found drug names: {drug_names}")
        return list(drug_names)
        
    except Exception as e:
        logger.error(f"Error reading drugnames.csv: {e}")
        return [drug_name]

def build_pubmed_query(drug_names: List[str]) -> str:
    """Build a PubMed query to search for case reports related to QT prolongation or TdP"""
    # Create drug name part of query with all variants from CSV
    drug_names_query = ' OR '.join(f'"{name}"[Title/Abstract]' for name in drug_names)
    
    # Build full query with cardiac terms
    query = (
        f'({drug_names_query}) AND '
        '(hERG[Title/Abstract] OR QT[Title/Abstract] OR QTc[Title/Abstract] OR '
        'torsad*[Title/Abstract]) AND '
        '"Humans"[Mesh] AND ("Case Reports"[Publication Type] OR "Cohort Studies"[Publication Type])'
    )
    
    logger.info(f"Searching PubMed: {query}")
    return query

def search_papers(drug_name: str, output_dir: Optional[str] = None) -> List[Dict]:
    """Search for papers about a drug and download PDFs."""
    # Set up output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'papers', drug_name.lower())
    os.makedirs(output_dir, exist_ok=True)
    
    # Get drug names from CSV
    drug_names = get_drug_names(drug_name)
    if not drug_names:
        logger.warning(f"No drug names found for {drug_name}")
        return []
        
    logger.info(f"Found drug names: {drug_names}")
    
    # Build and execute query
    query = build_pubmed_query(drug_names)
    logger.info(f"Searching PubMed: {query}")
    
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        logger.warning("No papers found")
        return []
    
    # Get paper details
    handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
    papers = list(Medline.parse(handle))
    handle.close()
    
    # Add PDF paths and extract text
    results = []
    for paper in papers:
        try:
            # Get DOI - check multiple fields
            doi = None
            # Check LID field
            if 'LID' in paper:
                for lid in paper['LID']:
                    if '[doi]' in lid:
                        doi = lid.replace('[doi]', '').strip()
                        break
            
            # Check AID field if no DOI found
            if not doi and 'AID' in paper:
                for aid in paper['AID']:
                    if '[doi]' in aid:
                        doi = aid.replace('[doi]', '').strip()
                        break
            
            # Check EID field if still no DOI
            if not doi and 'EID' in paper:
                doi = paper['EID'].strip()
            
            # Clean and validate DOI
            if doi:
                # Remove any URL prefixes
                doi = re.sub(r'^https?://doi.org/', '', doi)
                # Add DOI to paper
                paper['DOI'] = doi
                logger.info(f"Found DOI: {doi}")
            else:
                logger.warning(f"No DOI found for paper: {paper.get('TI', 'Unknown title')}")
                paper['DOI'] = ''
            
            # Try to get PDF
            pdf_filename = f"{paper.get('PMID', 'unknown')}.pdf"
            pdf_path = get_from_scihub(doi, output_dir, pdf_filename) if doi else None
            
            if pdf_path:
                # Convert PDF to text
                from .case_report_analyzer import convert_pdf_to_text
                text_path = convert_pdf_to_text(pdf_path)
                if text_path:
                    with open(text_path, 'r', encoding='utf-8') as f:
                        paper['full_text'] = f.read()
                
                paper['pdf_path'] = pdf_path
                results.append(paper)
            
        except Exception as e:
            logger.error(f"Error processing paper {paper.get('TI', 'Unknown title')}: {e}")
            continue
    
    logger.info(f"Found {len(results)} papers with PDFs")
    return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = search_papers(drug_name)
    print(f"\nFound {len(papers)} papers")

if __name__ == "__main__":
    main()
