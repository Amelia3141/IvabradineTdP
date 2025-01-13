import os
import sys
import time
import json
import logging
import pandas as pd
import requests
import io
import csv
import random
from bs4 import BeautifulSoup
import PyPDF2
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urljoin
from typing import List, Dict, Any, Optional, Tuple
from .case_report_analyzer import CaseReportAnalyzer
from Bio import Entrez

# Suppress SSL warnings
requests.packages.urllib3.disable_warnings()

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Set your email and API key for NCBI from environment variables
logger.info("Setting up NCBI credentials")
Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
Entrez.api_key = os.environ.get('NCBI_API_KEY')
logger.info(f"Using email: {Entrez.email}")
logger.info("NCBI credentials set up")

def get_full_text(pmid: str) -> Optional[str]:
    """Get full text for a paper using multiple methods."""
    try:
        # Try PMC first
        pmc_id = get_pmcid_from_pmid(pmid)
        if pmc_id:
            logger.info(f"Found PMC ID for {pmid}: {pmc_id}")
            text = get_pmc_text(pmc_id)
            if text and len(text) > 500:
                logger.info(f"Got {len(text)} chars from PMC for {pmid}")
                return text
        
        # Try Crossref next
        doi = get_doi_from_crossref(pmid)
        if doi:
            logger.info(f"Found DOI for {pmid}: {doi}")
            text = get_crossref_text(doi)
            if text and len(text) > 500:
                logger.info(f"Got {len(text)} chars from Crossref for {pmid}")
                return text
        
        # Try Sci-Hub as last resort
        logger.info(f"Trying Sci-Hub for {pmid}")
        text = get_scihub_text(pmid)
        if text and len(text) > 500:
            logger.info(f"Got {len(text)} chars from Sci-Hub for {pmid}")
            return text
            
        logger.warning(f"Could not get full text for {pmid} from any source")
        return None
        
    except Exception as e:
        logger.error(f"Error getting full text for {pmid}: {e}")
        return None

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

def get_pmc_text(pmcid: str) -> Optional[str]:
    """Get full text from PMC."""
    try:
        # Construct PMC URL
        url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        }
        
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code != 200:
            logger.warning(f"Got status code {response.status_code} from PMC")
            return None
            
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find the main article content
        article = soup.find('div', {'class': 'jig-ncbiinpagenav'})
        if not article:
            article = soup.find('article')
        if not article:
            article = soup.find('div', {'class': 'article'})
        
        if not article:
            logger.warning("Could not find article content in PMC")
            return None
            
        # Extract text, preserving some structure
        text_parts = []
        
        # Get title
        title = soup.find('h1', {'class': 'content-title'})
        if title:
            text_parts.append(title.get_text())
            
        # Get abstract
        abstract = soup.find('div', {'class': 'abstract'})
        if abstract:
            text_parts.append(abstract.get_text())
            
        # Get main content sections
        sections = article.find_all(['div', 'section'], {'class': ['section', 'sec']})
        for section in sections:
            # Get section title
            section_title = section.find(['h2', 'h3', 'h4'])
            if section_title:
                text_parts.append(section_title.get_text())
            
            # Get paragraphs
            paragraphs = section.find_all('p')
            for p in paragraphs:
                text_parts.append(p.get_text())
                
        text = '\n\n'.join(text_parts)
        text = ' '.join(text.split())  # Normalize whitespace
        
        if len(text) > 500:
            logger.info(f"Successfully extracted {len(text)} characters from PMC")
            return text
        else:
            logger.warning(f"PMC text too short: {len(text)} chars")
            return None
            
    except Exception as e:
        logger.error(f"Error getting PMC text: {e}")
        return None

def get_crossref_text(doi: str) -> Optional[str]:
    """Get full text using DOI from Crossref."""
    try:
        if not doi:
            return None
            
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5'
        }
        
        # Try to get the paper URL from Crossref
        crossref_url = f"https://api.crossref.org/works/{doi}"
        logger.info(f"Querying Crossref API: {crossref_url}")
        
        response = requests.get(crossref_url, headers=headers, timeout=30)
        if response.status_code != 200:
            logger.warning(f"Crossref API returned status code {response.status_code}")
            return None
            
        data = response.json()
        if not isinstance(data, dict) or 'message' not in data:
            logger.warning("Invalid response from Crossref API")
            return None
            
        message = data['message']
        if not isinstance(message, dict):
            logger.warning("Invalid message format in Crossref response")
            return None
            
        # Try to get the paper URL
        urls = []
        
        # Check for free PDF URL
        if 'link' in message and isinstance(message['link'], list):
            for link in message['link']:
                if isinstance(link, dict) and 'URL' in link:
                    content_type = link.get('content-type', '').lower()
                    if 'pdf' in content_type:
                        urls.append(link['URL'])
        
        # Check for alternative URLs
        if 'resource' in message and isinstance(message['resource'], dict):
            primary_url = message['resource'].get('primaryURL')
            if primary_url:
                urls.append(primary_url)
            
        # Add the main URL if available
        if 'URL' in message:
            urls.append(message['URL'])
            
        logger.info(f"Found {len(urls)} URLs from Crossref")
        
        # Try each URL
        for url in urls:
            try:
                logger.info(f"Trying URL: {url}")
                response = requests.get(url, headers=headers, timeout=30, verify=True)
                
                if response.status_code == 200:
                    content_type = response.headers.get('content-type', '').lower()
                    
                    if 'pdf' in content_type:
                        # Handle PDF
                        pdf_file = io.BytesIO(response.content)
                        pdf_reader = PyPDF2.PdfReader(pdf_file)
                        text = ""
                        for page in pdf_reader.pages:
                            text += page.extract_text()
                        text = text.strip()
                        if len(text) > 500:
                            logger.info(f"Successfully extracted {len(text)} characters from PDF")
                            return text
                    else:
                        # Handle HTML
                        soup = BeautifulSoup(response.text, 'html.parser')
                        
                        # Try different article content selectors
                        article_text = (
                            soup.find('div', {'class': 'fulltext-view'}) or
                            soup.find('article') or
                            soup.find('div', {'class': 'article-body'}) or
                            soup.find('div', {'class': 'main-content'}) or
                            soup.find('div', {'class': 'content'})
                        )
                        
                        if article_text:
                            text = article_text.get_text(separator=' ', strip=True)
                            text = ' '.join(text.split())
                            if len(text) > 500:
                                logger.info(f"Successfully extracted {len(text)} characters from HTML")
                                return text
                            
            except Exception as e:
                logger.error(f"Error processing URL {url}: {e}")
                continue
                
        logger.info(f"No text found in Crossref for {doi}")
        return None
        
    except Exception as e:
        logger.error(f"Error in get_crossref_text: {e}")
        return None

def get_scihub_text(pmid: str) -> Optional[str]:
    """Get full text from Sci-Hub."""
    try:
        # Only use .se domain
        domain = 'https://sci-hub.se'
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8'
        }
        
        # Add initial delay
        delay = random.uniform(5, 10)
        logger.info(f"Waiting {delay:.1f}s before Sci-Hub request")
        time.sleep(delay)
        
        url = f"{domain}/pubmed/{pmid}"
        logger.info(f"Trying Sci-Hub URL: {url}")
        
        response = requests.get(url, headers=headers, timeout=30, verify=False)
        if response.status_code != 200:
            logger.warning(f"Got status code {response.status_code} from Sci-Hub")
            return None
            
        soup = BeautifulSoup(response.text, 'html.parser')
        iframe = soup.find('iframe', {'id': 'pdf'})
        
        if not iframe:
            logger.warning("No PDF iframe found on Sci-Hub")
            return None
            
        pdf_url = iframe.get('src', '')
        if not pdf_url:
            logger.warning("No PDF URL found in iframe")
            return None
            
        if not pdf_url.startswith('http'):
            pdf_url = urljoin(domain, pdf_url)
        
        # Add delay before PDF request
        delay = random.uniform(5, 10)
        logger.info(f"Waiting {delay:.1f}s before PDF request")
        time.sleep(delay)
        
        logger.info(f"Downloading PDF from: {pdf_url}")
        pdf_response = requests.get(pdf_url, headers=headers, timeout=30, verify=False)
        
        if pdf_response.status_code == 200:
            try:
                pdf_file = io.BytesIO(pdf_response.content)
                pdf_reader = PyPDF2.PdfReader(pdf_file)
                text = ""
                for page in pdf_reader.pages:
                    text += page.extract_text()
                text = text.strip()
                if len(text) > 500:
                    logger.info(f"Got {len(text)} chars from PDF")
                    return text
                else:
                    logger.warning(f"PDF text too short: {len(text)} chars")
            except Exception as e:
                logger.error(f"Error extracting PDF text: {e}")
        else:
            logger.warning(f"Got status code {pdf_response.status_code} for PDF")
            
        return None
        
    except Exception as e:
        logger.error(f"Error in get_scihub_text: {e}")
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
            # Get PMID
            pmid = record['MedlineCitation']['PMID']
            
            # Get DOI
            doi = None
            if 'ArticleIdList' in record['PubmedData']:
                for id in record['PubmedData']['ArticleIdList']:
                    if id.attributes['IdType'] == 'doi':
                        doi = str(id)
                        break
            
            # Create result dict with key fields
            paper_dict = {
                'Title': record['MedlineCitation']['Article']['ArticleTitle'],
                'Abstract': record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                'PMID': pmid,
                'DOI': doi or '',
                'PubMed URL': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                'Year': record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'Age': '',
                'Sex': '',
                'Oral Dose (mg)': '',
                'QT (ms)': '',
                'QTc (ms)': '',
                'Heart Rate (bpm)': '',
                'Blood Pressure (mmHg)': '',
                'TdP': 'No'
            }
            results.append(paper_dict)
        except Exception as e:
            logger.error(f"Error processing record: {str(e)}")
            continue
            
    return pd.DataFrame(results)

def process_papers(pmids: List[str]) -> List[Dict]:
    """Process a list of papers in parallel."""
    try:
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = []
            for pmid in pmids:
                future = executor.submit(get_full_text, pmid)
                futures.append((pmid, future))
                
            texts = []
            for pmid, future in futures:
                try:
                    text = future.result(timeout=300)  # 5 minute timeout per paper
                    if text and len(text) > 500:
                        logger.info(f"Got {len(text)} chars for {pmid}")
                        texts.append({
                            'pmid': pmid,
                            'full_text': text
                        })
                    else:
                        logger.warning(f"No substantial text for {pmid}")
                        texts.append({
                            'pmid': pmid,
                            'full_text': None
                        })
                except Exception as e:
                    logger.error(f"Error processing {pmid}: {e}")
                    texts.append({
                        'pmid': pmid,
                        'full_text': None
                    })
            
        logger.info(f"Retrieved {len(texts)} full texts")
        return texts
        
    except Exception as e:
        logger.error(f"Error in process_papers: {e}")
        return []

def analyze_literature(drug_name: str) -> pd.DataFrame:
    """Analyze literature for a given drug."""
    try:
        # Search PubMed
        df = search_pubmed_case_reports(drug_name)
        if df.empty:
            logger.warning("No results found in PubMed")
            return pd.DataFrame()
            
        # Initialize all required columns
        required_columns = [
            'Title', 'Abstract', 'PMID', 'DOI', 'PubMed URL', 'Year',
            'Age', 'Sex', 'Oral Dose (mg)', 'QT (ms)', 'QTc (ms)',
            'Heart Rate (bpm)', 'Blood Pressure (mmHg)', 'TdP'
        ]
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = ''
        
        # Get PMIDs
        pmids = df['PMID'].tolist()
        logger.info(f"Found {len(pmids)} papers to analyze")
        
        # Process papers in parallel
        case_reports = process_papers(pmids)
        logger.info(f"Retrieved {len(case_reports)} full texts")
        
        # Analyze each case report
        for report in case_reports:
            try:
                pmid = report['pmid']
                text = report.get('full_text', '')
                
                if text and len(text) > 100:  # Ensure substantial text
                    logger.info(f"Text length for {pmid}: {len(text)} characters")
                    
                    analyzer = CaseReportAnalyzer()
                    analyzed = analyzer.analyze(text)
                    
                    # Update DataFrame with analyzed fields
                    idx = df.index[df['PMID'] == pmid].tolist()
                    if idx:
                        idx = idx[0]
                        df.loc[idx, 'Age'] = analyzed.get('age', '')
                        df.loc[idx, 'Sex'] = analyzed.get('sex', '')
                        df.loc[idx, 'Oral Dose (mg)'] = analyzed.get('oral_dose_value', '')
                        df.loc[idx, 'QT (ms)'] = analyzed.get('qt_value', '')
                        df.loc[idx, 'QTc (ms)'] = analyzed.get('qtc_value', '')
                        df.loc[idx, 'Heart Rate (bpm)'] = analyzed.get('heart_rate_value', '')
                        df.loc[idx, 'Blood Pressure (mmHg)'] = analyzed.get('blood_pressure_value', '')
                        df.loc[idx, 'TdP'] = 'Yes' if analyzed.get('tdp_present', False) else 'No'
                        
                        # Log what was found
                        found_fields = {k: v for k, v in analyzed.items() if v}
                        logger.info(f"Found fields for {pmid}: {found_fields}")
                    else:
                        logger.warning(f"Could not find PMID {pmid} in DataFrame")
                else:
                    logger.warning(f"No substantial text available for {pmid}")
                    
            except Exception as e:
                logger.error(f"Error analyzing paper {pmid}: {e}")
                continue
        
        # Sort by year descending
        if 'Year' in df.columns:
            df = df.sort_values('Year', ascending=False)
        
        # Reorder columns, keeping only those that exist
        columns = [col for col in required_columns if col in df.columns]
        df = df[columns]
        
        # Log summary of papers with data
        data_columns = ['Age', 'Sex', 'QT (ms)', 'QTc (ms)', 'Heart Rate (bpm)']
        data_columns = [col for col in data_columns if col in df.columns]
        if data_columns:
            found_data = df[df[data_columns].notna().any(axis=1)]
            logger.info(f"Found data in {len(found_data)} out of {len(df)} papers")
        
        return df
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {e}")
        return pd.DataFrame()

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers)} papers")

if __name__ == "__main__":
    main()
