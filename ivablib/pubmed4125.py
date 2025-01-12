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
try:
    import streamlit as st
    Entrez.email = st.secrets["NCBI_EMAIL"]
    Entrez.api_key = st.secrets["NCBI_API_KEY"]
except:
    Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
logger.info(f"Using email: {Entrez.email}")
logger.info("NCBI credentials set up")

def get_full_text(pmid, doi=None):
    """Get full text from PMC, Crossref, or Sci-Hub in that order. Falls back to abstract if full text unavailable."""
    text = None
    source = None
    logger.info(f"Attempting to get full text for PMID {pmid}, DOI {doi}")

    # Try PMC first
    try:
        pmc_id = get_pmcid_from_pmid(pmid)
        logger.info(f"PMC ID for {pmid}: {pmc_id}")
        if pmc_id:
            text = get_pmc_text(pmc_id)
            if text:
                source = 'PMC'
                logger.info(f"Got full text from PMC for {pmid}")
    except Exception as e:
        logger.warning(f"Error getting PMC text for {pmid}: {str(e)}")

    # Try Crossref if PMC failed
    if not text and doi:
        try:
            text = get_crossref_text(doi)
            if text:
                source = 'Crossref'
                logger.info(f"Got full text from Crossref for {pmid}")
        except Exception as e:
            logger.warning(f"Error getting Crossref text for {pmid}: {str(e)}")

    # Try Sci-Hub as last resort
    if not text:
        try:
            text = get_scihub_text(pmid, doi)
            if text:
                source = 'Sci-Hub'
                logger.info(f"Got full text from Sci-Hub for {pmid}")
        except Exception as e:
            logger.warning(f"Error getting Sci-Hub text for {pmid}: {str(e)}")

    # If still no text, try to get abstract from PubMed
    if not text:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
            abstract = handle.read()
            handle.close()
            if abstract:
                text = abstract
                source = 'PubMed Abstract'
                logger.info(f"Using PubMed abstract for {pmid} as fallback")
        except Exception as e:
            logger.warning(f"Error getting abstract for {pmid}: {str(e)}")

    if not text:
        logger.warning(f"No text returned for {pmid}")
        return None, None

    return text, source

def get_pmcid_from_pmid(pmid):
    """Convert PMID to PMCID using NCBI's ID converter."""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=tdp_analyzer&email={Entrez.email}&ids={pmid}&format=json"
        response = requests.get(url)
        data = response.json()
        if 'records' in data and data['records']:
            return data['records'][0].get('pmcid')
    except Exception as e:
        logger.error(f"Error converting PMID to PMCID: {e}")
    return None

def get_pmc_text(pmcid):
    """Get full text from PMC."""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        headers = {
            'User-Agent': 'Mozilla/5.0 TdPAnalyzer/1.0',
            'Accept': 'text/html,application/xhtml+xml'
        }
        response = requests.get(url, headers=headers)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Try multiple possible content containers in order
        article_text = None
        selectors = [
            ('div', {'class': 'jig-ncbiinpagenav'}),
            ('div', {'id': 'mc'}),  # Main content div in PMC
            ('div', {'class': 'article-body'}),
            ('div', {'id': 'body'}),
            ('article', {}),
            ('div', {'class': 'content'}),
            ('div', {'class': 'article'}),
        ]
        
        for tag, attrs in selectors:
            content = soup.find(tag, attrs)
            if content:
                # Remove references, supplementary material, etc.
                for div in content.find_all(['div', 'section'], {'class': ['ref-list', 'supplementary-material', 'copyright']}):
                    div.decompose()
                article_text = content
                break
        
        if article_text:
            # Clean the text
            text = article_text.get_text(separator=' ', strip=True)
            # Remove excessive whitespace
            text = ' '.join(text.split())
            return text
            
    except Exception as e:
        logger.error(f"Error getting PMC text: {e}")
    return None

def get_crossref_text(doi):
    """Get full text using DOI from Crossref."""
    try:
        # Try to get the paper URL from Crossref
        crossref_url = f"https://api.crossref.org/works/{doi}"
        response = requests.get(crossref_url)
        if response.status_code != 200:
            return None
            
        data = response.json()
        if 'message' not in data:
            return None
            
        # Get the paper URL
        paper_url = None
        if 'link' in data['message']:
            for link in data['message']['link']:
                if link.get('content-type', '').startswith('text/html'):
                    paper_url = link['URL']
                    break
                    
        if not paper_url and 'URL' in data['message']:
            paper_url = data['message']['URL']
            
        if not paper_url:
            return None
            
        # Get the paper content
        headers = {'User-Agent': 'Mozilla/5.0'}
        response = requests.get(paper_url, headers=headers)
        if response.status_code != 200:
            return None
            
        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Remove script and style elements
        for script in soup(["script", "style"]):
            script.decompose()
            
        # Get text and clean it up
        text = soup.get_text()
        
        # Break into lines and remove leading/trailing space
        lines = (line.strip() for line in text.splitlines())
        
        # Break multi-headlines into a line each
        chunks = (phrase.strip() for line in lines for phrase in line.split("  "))
        
        # Drop blank lines
        text = ' '.join(chunk for chunk in chunks if chunk)
        
        return text
        
    except Exception as e:
        logger.error(f"Error getting Crossref text: {str(e)}")
        return None

def get_scihub_text(pmid, doi=None):
    """Get full text from Sci-Hub."""
    try:
        # Try multiple Sci-Hub domains
        domains = [
            'https://sci-hub.se',
            'https://sci-hub.st',
            'https://sci-hub.ee',
            'https://sci-hub.wf',
            'https://sci-hub.ren'
        ]
        
        headers = {
            'User-Agent': 'Mozilla/5.0 TdPAnalyzer/1.0',
            'Accept': 'text/html,application/xhtml+xml'
        }
        
        for domain in domains:
            try:
                # Try DOI first if available
                if doi:
                    url = f"{domain}/{doi}"
                else:
                    url = f"{domain}/pubmed/{pmid}"
                
                response = requests.get(url, headers=headers, timeout=10)
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Look for embedded PDF
                pdf_iframe = soup.find('iframe', {'id': 'pdf'})
                if pdf_iframe and 'src' in pdf_iframe.attrs:
                    pdf_url = pdf_iframe['src']
                    if not pdf_url.startswith('http'):
                        pdf_url = f"https:{pdf_url}"
                    
                    # Download PDF
                    pdf_response = requests.get(pdf_url, headers=headers, timeout=10)
                    if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
                        pdf_file = io.BytesIO(pdf_response.content)
                        pdf_reader = PyPDF2.PdfReader(pdf_file)
                        text = ""
                        for page in pdf_reader.pages:
                            text += page.extract_text()
                        return text
            except Exception as e:
                continue
                
    except Exception as e:
        logger.error(f"Error getting Sci-Hub text: {e}")
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
    """Search PubMed for case reports about a drug"""
    try:
        # Get all related drug names
        drug_names = get_drug_names(drug_name)
        logger.info(f"Found {len(drug_names)} drug names for {drug_name}: {drug_names}")
        
        # Build search query
        drug_query = ' OR '.join([f'"{name}"[Title/Abstract]' for name in drug_names])
        query = f'({drug_query}) AND (hERG[Title/Abstract] OR QT[Title/Abstract] OR QTc[Title/Abstract] OR torsad*[Title/Abstract]) AND "Humans"[Mesh] AND ("Case Reports"[Publication Type])'
        
        logger.info(f"Searching PubMed: {query}")
        
        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
            logger.info("Got esearch handle")
            record = Entrez.read(handle)
            logger.info(f"Search returned {len(record['IdList'])} results")
            handle.close()
            
            if not record["IdList"]:
                logger.info("No results found")
                return pd.DataFrame()
                
            # Fetch details
            logger.info(f"Fetching details for {len(record['IdList'])} papers")
            handle = Entrez.efetch(db="pubmed", id=record["IdList"], 
                                 rettype="xml", retmode="text")
            records = Entrez.read(handle)
            logger.info(f"Got {len(records['PubmedArticle'])} paper details")
            handle.close()

            results = []
            for paper in records['PubmedArticle']:
                try:
                    paper_dict = {
                        'PMID': paper['MedlineCitation']['PMID'],
                        'Title': paper['MedlineCitation']['Article']['ArticleTitle'],
                        'Abstract': paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                        'Year': paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', '')
                    }
                    results.append(paper_dict)
                except Exception as e:
                    logger.error(f"Error processing paper: {e}")
                    continue
            
            logger.info(f"Processed {len(results)} papers successfully")
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            return pd.DataFrame()
            
    except Exception as e:
        logger.error(f"Error in search_pubmed: {e}")
        return pd.DataFrame()

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
                                     rettype="xml", retmode="text")
                records = Entrez.read(handle)
                handle.close()

                return records['PubmedArticle']
                
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
            paper_dict = {
                'PMID': record['MedlineCitation']['PMID'],
                'Case Report Title': record['MedlineCitation']['Article']['ArticleTitle'],
                'Abstract': record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                'Year': record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'Age': '',  # Will be extracted by regex
                'Sex': '',  # Will be extracted by regex
                'Oral Dose (mg)': '',  # Will be extracted by regex
                'Theoretical Max Concentration (μM)': '',  # Will be extracted by regex
                '40% Bioavailability': '',  # Will be extracted by regex
                'Theoretical HERG IC50 / Concentration μM': '',  # Will be extracted by regex
                '40% Plasma Concentration': '',  # Will be extracted by regex
                'Uncorrected QT (ms)': '',  # Will be extracted by regex
                'QTc': '',  # Will be extracted by regex
                'QTR': '',  # Will be calculated later
                'QTF': '',  # Will be calculated later
                'Heart Rate (bpm)': '',  # Will be extracted by regex
                'Torsades de Pointes?': 'No',  # Will be updated by regex
                'Blood Pressure (mmHg)': '',  # Will be extracted by regex
                'Medical History': '',  # Will be extracted by regex
                'Medication History': '',  # Will be extracted by regex
                'Course of Treatment': ''  # Will be extracted by regex
            }
            results.append(paper_dict)
        except Exception as e:
            logger.error(f"Error processing record: {str(e)}")
            continue
            
    return pd.DataFrame(results)

def get_texts_parallel(pmids):
    """Get full texts for multiple PMIDs in parallel using ThreadPoolExecutor"""
    texts = {}
    logger.info(f"Getting full texts for {len(pmids)} papers")
    
    # First get DOIs for all papers
    dois = {}
    for pmid in pmids:
        try:
            # Query NCBI for DOI
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
            records = Entrez.read(handle)
            handle.close()
            
            article = records['PubmedArticle'][0]['MedlineCitation']['Article']
            # Try both ELocationID and ArticleId for DOI
            doi = None
            
            # Check ELocationID
            if 'ELocationID' in article:
                for id in article['ELocationID']:
                    if id.attributes['EIdType'] == 'doi':
                        doi = str(id)
                        break
            
            # If not found, check ArticleIdList
            if not doi and 'ArticleIdList' in records['PubmedArticle'][0]['PubmedData']:
                for id in records['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
                    if id.attributes['IdType'] == 'doi':
                        doi = str(id)
                        break
            
            if doi:
                dois[pmid] = doi
                logger.info(f"Found DOI for {pmid}: {doi}")
            
            # Add a small delay between requests
            time.sleep(1)
            
        except Exception as e:
            logger.error(f"Error getting DOI for {pmid}: {e}")
    
    def get_text_with_delay(pmid):
        """Wrapper to add delay between requests"""
        text, source = get_full_text(pmid, dois.get(pmid))
        time.sleep(2)  # Add delay between full text requests
        return text, source
    
    # Use ThreadPoolExecutor to run requests in parallel
    with ThreadPoolExecutor(max_workers=3) as executor:  # Reduced max_workers to avoid overwhelming servers
        future_to_pmid = {
            executor.submit(get_text_with_delay, pmid): pmid 
            for pmid in pmids
        }
        
        for future in concurrent.futures.as_completed(future_to_pmid):
            pmid = future_to_pmid[future]
            try:
                result = future.result()
                if result:
                    text, source = result
                    if text:
                        texts[pmid] = text
                        logger.info(f"Successfully got text for {pmid} from {source}")
                    else:
                        logger.warning(f"No text returned for {pmid}")
            except Exception as e:
                logger.error(f"Error processing {pmid}: {e}")
                
    logger.info(f"Got {len(texts)} full texts out of {len(pmids)} papers")
    return texts

def analyze_literature(drug_name):
    """Analyze literature for QT-related case reports."""
    try:
        logger.info(f"Starting literature analysis for {drug_name}")
        
        # Search for papers
        papers = search_pubmed_case_reports(drug_name)
        if not papers:
            logger.info(f"No papers found for {drug_name}")
            return {
                'case_reports': [],
                'message': f"No case reports found for {drug_name}"
            }
            
        # Get full texts
        papers = get_full_texts(papers)
        if not papers:
            logger.info(f"No full texts found for {drug_name}")
            return {
                'case_reports': [],
                'message': f"No full texts available for {drug_name}"
            }
        
        # Process papers
        result = process_papers(papers)
        if 'error' in result:
            return result
        
        if not result.get('case_reports'):
            return {
                'case_reports': [],
                'message': f"No relevant case reports found for {drug_name}"
            }
            
        return result
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {str(e)}")
        return {
            'error': f"Error analyzing literature: {str(e)}"
        }

def get_full_texts(papers):
    """Get full texts for papers"""
    papers_with_text = []
    total_papers = len(papers)
    texts_found = 0
    
    for paper in papers:
        pmid = str(paper.get('PMID', ''))
        doi = paper.get('DOI', '')
        
        logger.info(f"Getting text for paper {pmid}")
        
        # Get text from various sources
        text, source = get_full_text(pmid, doi)
        
        if text:
            texts_found += 1
            paper_with_text = {
                'PMID': pmid,
                'Title': paper.get('Title', ''),
                'Abstract': paper.get('Abstract', ''),
                'FullText': text,
                'TextSource': source,
                'DOI': doi
            }
            papers_with_text.append(paper_with_text)
            logger.info(f"Added text for paper {pmid} from {source}")
        else:
            logger.warning(f"No text found for paper {pmid}")
    
    logger.info(f"Found text for {texts_found} out of {total_papers} papers")
    return papers_with_text

def process_papers(papers):
    """Process papers to extract relevant information"""
    try:
        analyzer = CaseReportAnalyzer()
        results = []
        
        for paper in papers:
            # Extract information using the analyzer
            info = analyzer.analyze_paper(paper)
            if info:
                results.append(info)
                logger.info(f"Successfully processed paper {paper.get('PMID', 'Unknown')}")
            else:
                logger.warning(f"No relevant information found in paper {paper.get('PMID', 'Unknown')}")
                
        if results:
            logger.info(f"Successfully processed {len(results)} papers")
            return {'case_reports': results}
        else:
            logger.warning("No case reports found in any papers")
            return {'message': "No case reports found in the literature."}
        
    except Exception as e:
        logger.error(f"Error processing papers: {str(e)}")
        return {'error': str(e)}

def get_paper_details(pmids):
    """Get paper details from PubMed"""
    try:
        # Fetch details
        logger.info(f"Fetching details for {len(pmids)} papers")
        handle = Entrez.efetch(db="pubmed", id=pmids, 
                             rettype="xml", retmode="text")
        records = Entrez.read(handle)
        logger.info(f"Got {len(records['PubmedArticle'])} paper details")
        handle.close()

        results = []
        for paper in records['PubmedArticle']:
            try:
                paper_dict = {
                    'PMID': paper['MedlineCitation']['PMID'],
                    'Title': paper['MedlineCitation']['Article']['ArticleTitle'],
                    'Abstract': paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                    'Year': paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', '')
                }
                results.append(paper_dict)
            except Exception as e:
                logger.error(f"Error processing paper: {e}")
                continue
        
        logger.info(f"Processed {len(results)} papers successfully")
        return pd.DataFrame(results)
        
    except Exception as e:
        logger.error(f"Error getting paper details: {e}")
        return pd.DataFrame()

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers['case_reports'])} papers")

if __name__ == "__main__":
    main()
