import os
import logging
import pandas as pd
import requests
from bs4 import BeautifulSoup
import re
import json
from Bio import Entrez
import PyPDF2
import io
import logging
from urllib.parse import urljoin
from typing import List, Dict, Optional
import sys
import time
import difflib
import csv
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

# Suppress SSL warnings
requests.packages.urllib3.disable_warnings()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email and API key for NCBI from environment variables
logger.info("Setting up NCBI credentials")
Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
Entrez.api_key = os.environ.get('NCBI_API_KEY')
logger.info(f"Using email: {Entrez.email}")
logger.info("NCBI credentials set up")

def get_full_text(pmid, doi=None):
    """Get full text from PMC, Crossref, or Sci-Hub in that order."""
    text = None
    source = None

    # Try PMC first
    try:
        pmc_id = get_pmcid_from_pmid(pmid)
        if pmc_id:
            text = get_pmc_text(pmc_id)
            if text:
                source = 'PMC'
                return text, source
    except Exception as e:
        logging.error(f"Error getting PMC text for {pmid}: {e}")

    # Try Crossref next if we have a DOI
    if doi:
        try:
            crossref_text = get_crossref_text(doi)
            if crossref_text:
                text = crossref_text
                source = 'Crossref'
                return text, source
        except Exception as e:
            logging.error(f"Error getting Crossref text for {doi}: {e}")

    # Finally try Sci-Hub
    try:
        scihub_text = get_scihub_text(pmid, doi)
        if scihub_text:
            text = scihub_text
            source = 'Sci-Hub'
            return text, source
    except Exception as e:
        logging.error(f"Error getting Sci-Hub text for {pmid}: {e}")

    return text, source

def get_pmcid_from_pmid(pmid):
    """Convert PMID to PMCID using NCBI's ID converter."""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=my.email@example.com&ids={pmid}&format=json"
        response = requests.get(url)
        data = response.json()
        if 'records' in data and data['records']:
            return data['records'][0].get('pmcid')
    except Exception as e:
        logging.error(f"Error converting PMID to PMCID: {e}")
    return None

def get_pmc_text(pmcid):
    """Get full text from PMC."""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        article_text = soup.find('div', {'class': 'jig-ncbiinpagenav'})
        if article_text:
            return article_text.get_text()
    except Exception as e:
        logging.error(f"Error getting PMC text: {e}")
    return None

def get_crossref_text(doi):
    """Get full text using DOI from Crossref."""
    try:
        url = f"https://api.crossref.org/works/{doi}"
        response = requests.get(url)
        data = response.json()
        if 'message' in data and 'URL' in data['message']:
            publisher_url = data['message']['URL']
            response = requests.get(publisher_url)
            soup = BeautifulSoup(response.text, 'html.parser')
            article_text = soup.find('div', {'class': ['article-body', 'main-content']})
            if article_text:
                return article_text.get_text()
    except Exception as e:
        logging.error(f"Error getting Crossref text: {e}")
    return None

def get_scihub_text(pmid, doi=None):
    """Get full text from Sci-Hub."""
    try:
        identifier = doi if doi else pmid
        url = f"https://sci-hub.se/{identifier}"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        pdf_url = soup.find('iframe', {'id': 'pdf'})
        if pdf_url:
            pdf_url = urljoin('https://sci-hub.se/', pdf_url['src'])
            pdf_response = requests.get(pdf_url)
            pdf_file = io.BytesIO(pdf_response.content)
            pdf_reader = PyPDF2.PdfReader(pdf_file)
            text = ""
            for page in pdf_reader.pages:
                text += page.extract_text()
            return text
    except Exception as e:
        logging.error(f"Error getting Sci-Hub text: {e}")
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
        '"Humans"[Mesh] AND ("Case Reports"[Publication Type])'
    )
    
    logger.info(f"Searching PubMed: {query}")
    return query

def search_pubmed_case_reports(drug_name: str) -> pd.DataFrame:
    """Search PubMed for case reports about the drug."""
    try:
        query = build_pubmed_query(get_drug_names(drug_name))
        records = _search_pubmed(query)
        return format_pubmed_results(records)
    except Exception as e:
        logger.error(f"Error in case reports search: {str(e)}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

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
                records = list(Entrez.parse(handle, 'medline'))
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

def get_texts_parallel(pmids):
    """Get full texts for multiple PMIDs in parallel using ThreadPoolExecutor"""
    texts = {}
    
    def fetch_single_text(pmid):
        try:
            # Get DOI from PubMed
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
            xml = handle.read()
            handle.close()
            
            doi_match = re.search(r'<ELocationID EIdType="doi">(.*?)</ELocationID>', str(xml))
            doi = doi_match.group(1) if doi_match else None
            
            text, source = get_full_text(pmid, doi)
            if text:
                logger.info(f"Got text for {pmid} from {source}")
                return pmid, text
            return None
        except Exception as e:
            logger.error(f"Error fetching text for {pmid}: {e}")
            return None

    # Use ThreadPoolExecutor to run requests in parallel
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(fetch_single_text, pmid) for pmid in pmids]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                pmid, text = result
                texts[pmid] = text
                
    return texts

def analyze_literature(drug_name: str) -> Dict:
    """Analyze literature for a drug"""
    try:
        logger.info(f"Searching literature for {drug_name}")
        
        # Only search for case reports
        case_reports = search_pubmed_case_reports(drug_name)
        
        if not case_reports.empty:
            # Get PMIDs
            pmids = [p for p in case_reports['PMID'].tolist() if p != 'No PMID']
            
            # Get texts in parallel using ThreadPoolExecutor
            texts = get_texts_parallel(pmids)
            
            # Create paper objects
            for _, paper in case_reports.iterrows():
                pmid = paper['PMID']
                if pmid in texts:
                    papers_to_analyze.append({
                        'pmid': pmid,
                        'title': paper['Title'],
                        'authors': paper['Authors'].split(', '),
                        'year': paper['Year'],
                        'journal': paper['Journal'],
                        'full_text': texts[pmid]
                    })
        
        # Analyze papers using CaseReportAnalyzer
        from .case_report_analyzer import analyze_papers
        analysis_results = analyze_papers(papers_to_analyze, drug_name)
        
        return {
            "case_reports": case_reports.to_dict('records'),
            "full_texts": papers_to_analyze,
            "analysis": analysis_results.to_dict('records') if not analysis_results.empty else []
        }
        
    except Exception as e:
        logger.error(f"Error analyzing literature: {str(e)}")
        return {"error": str(e)}

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers['case_reports'])} papers")

if __name__ == "__main__":
    main()
