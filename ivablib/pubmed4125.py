import os
import logging
import pandas as pd
import requests
import urllib3
from bs4 import BeautifulSoup
try:
    from Bio import Entrez
    from Bio import Medline
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False
    logging.warning("Biopython not found. Literature analysis will be limited.")

from typing import List, Dict, Optional
import sys
import time
import re
from urllib.parse import quote
try:
    import PyPDF2
    HAVE_PYPDF2 = True
except ImportError:
    HAVE_PYPDF2 = False
    logging.warning("PyPDF2 not found. PDF analysis will be limited.")

import difflib
import csv

# Suppress SSL warnings
urllib3.disable_warnings()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email and API key for NCBI from environment variables
logger.info("Setting up NCBI credentials")
if HAVE_BIOPYTHON:
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    logger.info(f"Using email: {Entrez.email}")
    logger.info("NCBI credentials set up")

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
                response = requests.get(sci_hub_url, verify=False, timeout=30)
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
                    pdf_response = requests.get(pdf_url, verify=False, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
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
                    pdf_response = requests.get(pdf_url, verify=False, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
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
        drug_names = set()
        drug_name_lower = drug_name.lower()
        csv_path = os.path.join(os.path.dirname(__file__), 'drugnames.csv')
        
        with open(csv_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            if not header or ',' not in header:
                header = "DrugName,ActiveIngredient"
            
            reader = csv.DictReader(f, fieldnames=['DrugName', 'ActiveIngredient'])
            
            for row in reader:
                if not row['DrugName'] or not row['ActiveIngredient']:
                    continue
                    
                brand_name = row['DrugName'].strip().lower()
                active_ingredient = row['ActiveIngredient'].strip().lower()
                
                # Check if the drug name matches either brand name or active ingredient
                if (drug_name_lower in brand_name or 
                    drug_name_lower in active_ingredient or
                    brand_name in drug_name_lower or
                    active_ingredient in drug_name_lower):
                    
                    # Skip combination products for better search focus
                    if 'and' in brand_name.lower() or ';' in active_ingredient.lower():
                        continue
                        
                    drug_names.add(row['DrugName'].strip())
                    drug_names.add(row['ActiveIngredient'].strip())
        
        # Remove empty strings and duplicates
        all_names = set()
        for name in drug_names:
            if name and ';' not in name and ' and ' not in name.lower():
                all_names.add(name)
        
        # Also add the original search term
        all_names.add(drug_name)
        
        # Remove empty strings and duplicates
        result = sorted(list(name for name in all_names if name))
        logger.info(f"Found {len(result)} drug names for {drug_name}: {result}")
        return result
        
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

def _search_pubmed(query: str) -> List[Dict]:
    """Base function to search PubMed"""
    try:
        # Set up retries
        max_attempts = 3
        for attempt in range(max_attempts):
            try:
                # Search PubMed
                handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
                record = Entrez.read(handle)
                handle.close()

                if not record["IdList"]:
                    return []

                # Fetch details
                handle = Entrez.efetch(db="pubmed", id=record["IdList"], 
                                     rettype="medline", retmode="text")
                records = list(Medline.parse(handle))
                handle.close()

                return records

            except Exception as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(2 ** attempt)  # Exponential backoff
                continue

    except Exception as e:
        logger.error(f"Error searching PubMed: {str(e)}")
        return []

def search_pubmed_case_reports(drug_name: str) -> pd.DataFrame:
    """Search PubMed for case reports about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND "case reports"[Publication Type] AND ("torsade de pointes" OR "qt prolongation")'
        records = _search_pubmed(query)
        return format_pubmed_results(records)
    except Exception as e:
        logger.error(f"Error in case reports search: {str(e)}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def search_pubmed_cohort_studies(drug_name: str) -> pd.DataFrame:
    """Search PubMed for cohort studies about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND ("cohort studies"[MeSH Terms] OR "cohort"[Title/Abstract]) AND ("torsade de pointes" OR "qt prolongation")'
        records = _search_pubmed(query)
        return format_pubmed_results(records)
    except Exception as e:
        logger.error(f"Error in cohort studies search: {str(e)}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def search_pubmed_clinical_trials(drug_name: str) -> pd.DataFrame:
    """Search PubMed for clinical trials about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND "clinical trial"[Publication Type] AND ("torsade de pointes" OR "qt prolongation")'
        records = _search_pubmed(query)
        return format_pubmed_results(records)
    except Exception as e:
        logger.error(f"Error in clinical trials search: {str(e)}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def format_pubmed_results(records: List[Dict]) -> pd.DataFrame:
    """Format PubMed results into a DataFrame."""
    results = []
    for record in records:
        try:
            # Extract authors
            authors = record.get("AU", [])
            author_str = ", ".join(authors) if authors else "No authors listed"
            
            # Get year from date
            date = record.get("DP", "")
            year = date.split()[0] if date else "N/A"
            
            results.append({
                "Title": record.get("TI", "No title available"),
                "Authors": author_str,
                "Journal": record.get("JT", "Journal not specified"),
                "Year": year,
                "PMID": record.get("PMID", "No PMID")
            })
        except Exception as e:
            logger.error(f"Error processing record: {str(e)}")
            continue
            
    return pd.DataFrame(results)

def analyze_literature(drug_name: str) -> Dict:
    """Analyze literature for a drug"""
    if not HAVE_BIOPYTHON:
        return {
            "error": "Biopython not installed. Please install with: pip install biopython"
        }
        
    try:
        logger.info(f"Searching literature for {drug_name}")
        
        # Get all related drug names
        drug_names = get_drug_names(drug_name)
        logger.info(f"Found {len(drug_names)} drug names for {drug_name}: {drug_names}")
        
        # Build search query with all drug names
        query = build_pubmed_query(drug_names)
        logger.info(f"Searching PubMed: {query}")
        
        # Search different types of studies
        case_reports = search_pubmed_case_reports(drug_name)
        cohort_studies = search_pubmed_cohort_studies(drug_name)
        clinical_trials = search_pubmed_clinical_trials(drug_name)
        
        # Get full text papers and analyze them
        papers_to_analyze = []
        for paper_type in [case_reports, cohort_studies, clinical_trials]:
            if not paper_type.empty:
                for _, paper in paper_type.iterrows():
                    try:
                        pmid = paper['PMID']
                        if pmid == 'No PMID':
                            continue
                            
                        # Get PDF from Sci-Hub
                        pdf_path = get_from_scihub(pmid, os.getcwd(), f"{pmid}.pdf")
                        if pdf_path:
                            text = extract_text_from_pdf(pdf_path)
                            if text:
                                papers_to_analyze.append({
                                    'pmid': pmid,
                                    'title': paper['Title'],
                                    'authors': paper['Authors'].split(', '),
                                    'year': paper['Year'],
                                    'journal': paper['Journal'],
                                    'full_text': text
                                })
                    except Exception as e:
                        logger.error(f"Error getting full text for {pmid}: {str(e)}")
                        continue
        
        # Analyze papers using CaseReportAnalyzer
        from .case_report_analyzer import analyze_papers
        analysis_results = analyze_papers(papers_to_analyze, drug_name)
        
        return {
            "case_reports": case_reports.to_dict('records'),
            "cohort_studies": cohort_studies.to_dict('records'),
            "clinical_trials": clinical_trials.to_dict('records'),
            "full_texts": papers_to_analyze,
            "analysis": analysis_results.to_dict('records') if not analysis_results.empty else []
        }
        
    except Exception as e:
        logger.error(f"Error analyzing literature: {str(e)}")
        return {"error": str(e)}

def search_papers(query: str, email: str, api_key: Optional[str] = None) -> List[Dict]:
    """Search PubMed for papers matching query"""
    try:
        if HAVE_BIOPYTHON:
            Entrez.email = email
            if api_key:
                Entrez.api_key = api_key
        
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()

            if not record['IdList']:
                return []

            # Fetch details for each paper
            papers = []
            for pmid in record['IdList']:
                try:
                    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                    record = Medline.read(handle)
                    handle.close()
                    
                    if not isinstance(record, dict):
                        continue
                        
                    paper = {
                        'pmid': pmid,
                        'title': record.get('TI', ''),
                        'authors': record.get('AU', []),
                        'journal': record.get('JT', ''),
                        'year': record.get('DP', '').split()[0] if record.get('DP') else '',
                        'abstract': record.get('AB', '')
                    }
                    papers.append(paper)
                    
                except Exception as e:
                    logger.error(f"Error fetching paper {pmid}: {str(e)}")
                    continue
                    
        else:
            logger.error("Biopython not installed. Please install with: pip install biopython")
            return []
        
        return papers
        
    except Exception as e:
        logger.error(f"Error searching PubMed: {str(e)}")
        return []

def search_papers(drug_name: str, output_dir: Optional[str] = None) -> List[Dict]:
    """Search for papers about a drug and download PDFs.
    
    Args:
        drug_name: Name of drug to search for
        output_dir: Directory to save PDFs to. If None, will use default.
        
    Returns:
        List of papers with their metadata and PDF paths
    """
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
    
    max_retries = 3
    retry_delay = 1  # seconds
    
    for attempt in range(max_retries):
        try:
            papers = search_papers(query, Entrez.email, Entrez.api_key)
            break
        except Exception as e:
            if attempt == max_retries - 1:
                logger.error(f"Failed to search PubMed after {max_retries} attempts: {e}")
                return []
            time.sleep(retry_delay)
            retry_delay *= 2  # exponential backoff
    
    if not papers:
        logger.warning("No papers found")
        return []
    
    # Add PDF paths
    results = []
    for paper in papers:
        try:
            # Get DOI
            doi = None
            if 'LID' in paper:
                doi = paper['LID'][0] if paper['LID'][0].endswith('[doi]') else None
            if not doi and 'AID' in paper:
                for aid in paper['AID']:
                    if aid.endswith('[doi]'):
                        doi = aid
                        break
            
            if not doi:
                logger.warning(f"No DOI found for paper: {paper.get('TI', 'Unknown title')}")
                continue
                
            # Clean DOI
            doi = doi.replace('[doi]', '').strip()
            
            # Try to get PDF
            pdf_filename = f"{paper.get('PMID', 'unknown')}.pdf"
            pdf_path = get_from_scihub(doi, output_dir, pdf_filename)
            
            if pdf_path:
                paper['pdf_path'] = pdf_path
                results.append(paper)
            
        except Exception as e:
            logger.error(f"Error processing paper {paper.get('TI', 'Unknown title')}: {e}")
            continue
    
    logger.info(f"Found {len(results)} papers with PDFs")
    return results

def analyze_text_for_findings(text: str) -> List[str]:
    """Analyze text for key findings related to TdP risk"""
    findings = []
    
    # Look for mentions of QT intervals
    qt_pattern = r'QTc?\s*(?:interval)?\s*(?:was|increased|prolonged|of)\s*(\d+(?:\.\d+)?)\s*(?:ms|msec)'
    qt_matches = re.finditer(qt_pattern, text, re.IGNORECASE)
    for match in qt_matches:
        qt_value = float(match.group(1))
        if qt_value > 440:  # Clinically significant threshold
            findings.append(f"QTc prolongation found: {qt_value} ms")
            
    # Look for TdP events
    tdp_pattern = r'(?:torsade\s+de\s+pointes|TdP|torsades)\s+(?:was|were|occurred|developed|observed)'
    if re.search(tdp_pattern, text, re.IGNORECASE):
        findings.append("TdP event reported")
        
    # Look for drug concentrations
    conc_pattern = r'(?:concentration|level)\s*(?:of|was|were)\s*(\d+(?:\.\d+)?)\s*(ng/ml|µg/ml|mg/l)'
    conc_matches = re.finditer(conc_pattern, text, re.IGNORECASE)
    for match in conc_matches:
        findings.append(f"Drug concentration: {match.group(1)} {match.group(2)}")
        
    return findings

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = search_papers(drug_name)
    print(f"\nFound {len(papers)} papers")

if __name__ == "__main__":
    main()
