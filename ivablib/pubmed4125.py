import os
import logging
import pandas as pd
import requests
from bs4 import BeautifulSoup
from Bio import Entrez
import time
import re
import io
import PyPDF2
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urljoin
from typing import List, Dict, Any, Optional, Tuple
from .case_report_analyzer import CaseReportAnalyzer
import random
import sys

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
    logger.info(f"Attempting to get full text for PMID {pmid}, DOI {doi}")

    # Try PMC first
    try:
        pmc_id = get_pmcid_from_pmid(pmid)
        logger.info(f"PMC ID for {pmid}: {pmc_id}")
        if pmc_id:
            text = get_pmc_text(pmc_id)
            if text:
                source = 'PMC'
                logger.info(f"Successfully got text from PMC for {pmid}")
                return text, source
            else:
                logger.info(f"No text found in PMC for {pmid}")
    except Exception as e:
        logger.error(f"Error getting PMC text for {pmid}: {e}")

    # Try Crossref next if we have a DOI
    if doi:
        try:
            logger.info(f"Trying Crossref for {pmid} with DOI {doi}")
            crossref_text = get_crossref_text(doi)
            if crossref_text:
                text = crossref_text
                source = 'Crossref'
                logger.info(f"Successfully got text from Crossref for {pmid}")
                return text, source
            else:
                logger.info(f"No text found in Crossref for {pmid}")
        except Exception as e:
            logger.error(f"Error getting Crossref text for {doi}: {e}")

    # Finally try Sci-Hub
    try:
        logger.info(f"Trying Sci-Hub for {pmid}")
        scihub_text = get_scihub_text(pmid, doi)
        if scihub_text:
            text = scihub_text
            source = 'Sci-Hub'
            logger.info(f"Successfully got text from Sci-Hub for {pmid}")
            return text, source
        else:
            logger.info(f"No text found in Sci-Hub for {pmid}")
    except Exception as e:
        logger.error(f"Error getting Sci-Hub text for {pmid}: {e}")

    logger.warning(f"Could not get full text for {pmid} from any source")
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
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
        }
        response = requests.get(url, headers=headers)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Try multiple possible content containers in order
        article_text = None
        
        # First try the main article content
        content = soup.find('article')
        if content:
            # Only remove specific sections that aren't part of the main text
            for div in content.find_all(['div', 'section'], {'class': ['ref-list', 'copyright']}):
                div.decompose()
            article_text = content
        
        # If no article tag, try other containers
        if not article_text:
            selectors = [
                ('div', {'class': 'jig-ncbiinpagenav'}),
                ('div', {'id': 'mc'}),
                ('div', {'class': 'article-body'}),
                ('div', {'id': 'body'}),
                ('div', {'class': 'content'}),
                ('div', {'class': 'article'})
            ]
            
            for tag, attrs in selectors:
                content = soup.find(tag, attrs)
                if content:
                    # Only remove specific sections
                    for div in content.find_all(['div', 'section'], {'class': ['ref-list', 'copyright']}):
                        div.decompose()
                    article_text = content
                    break
        
        if article_text:
            # Extract text while preserving some structure
            paragraphs = []
            for p in article_text.find_all(['p', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6']):
                text = p.get_text(strip=True)
                if text:  # Only add non-empty paragraphs
                    paragraphs.append(text)
            
            # Join paragraphs with newlines to preserve structure
            text = '\n'.join(paragraphs)
            
            # Clean up the text
            text = ' '.join(text.split())  # Remove excessive whitespace
            text = text.replace('- ', '')  # Remove hyphenation
            
            if len(text) > 500:  # Only return if we got substantial text
                logger.info(f"Successfully extracted {len(text)} characters from PMC")
                return text
            else:
                logger.warning("PMC text too short, might be incomplete")
            
        logger.warning("No article text found in PMC")
        
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
        if 'message' not in data:
            logger.warning("No message in Crossref response")
            return None
            
        # Try to get the paper URL
        urls = []
        
        # Check for free PDF URL
        if 'link' in data['message']:
            for link in data['message']['link']:
                if link.get('content-type', '').lower() == 'application/pdf':
                    urls.append(link['URL'])
        
        # Check for alternative URLs
        if 'resource' in data['message']:
            urls.extend([r['URL'] for r in data['message']['resource']])
            
        # Add the main URL if available
        if 'url' in data['message']:
            urls.append(data['message']['url'])
            
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
                
        logger.warning(f"Could not extract text from any Crossref URLs for DOI {doi}")
        return None
        
    except Exception as e:
        logger.error(f"Error in get_crossref_text: {e}")
        return None

def get_scihub_text(pmid, doi=None):
    """Get full text from Sci-Hub."""
    try:
        # Try multiple Sci-Hub domains
        domains = [
            'https://sci-hub.ru',
            'https://sci-hub.st',
            'https://sci-hub.se',
            'https://sci-hub.wf',
            'https://sci-hub.cat'
        ]
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        }
        
        session = requests.Session()
        session.headers.update(headers)
        session.verify = False  # Disable SSL verification
        
        # Add delay between requests
        time.sleep(random.uniform(2, 5))
        
        for domain in domains:
            try:
                logger.info(f"Trying Sci-Hub domain: {domain}")
                
                # Try DOI first if available, then PMID
                if doi:
                    url = f"{domain}/{doi}"
                else:
                    url = f"{domain}/pubmed/{pmid}"
                
                # Add delay between attempts
                time.sleep(random.uniform(1, 3))
                
                response = session.get(url, timeout=30)
                if response.status_code != 200:
                    logger.warning(f"Got status code {response.status_code} from {domain}")
                    continue
                
                if 'location not found' in response.text.lower() or 'article not found' in response.text.lower():
                    logger.info(f"Paper not found on {domain}")
                    continue
                    
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # First try to get text directly from the page
                article_text = soup.find('div', {'id': 'pdf'}) or soup.find('embed', {'id': 'pdf'})
                if article_text:
                    text = article_text.get_text(separator=' ', strip=True)
                    if len(text) > 500:
                        logger.info(f"Got {len(text)} chars directly from page")
                        return text
                
                # Look for embedded PDF
                pdf_iframe = soup.find('iframe', {'id': 'pdf'}) or soup.find('embed', {'id': 'pdf'})
                if pdf_iframe and 'src' in pdf_iframe.attrs:
                    pdf_url = pdf_iframe['src']
                    if not pdf_url.startswith('http'):
                        pdf_url = urljoin(domain, pdf_url) if pdf_url.startswith('/') else f"https:{pdf_url}"
                    
                    logger.info(f"Found PDF at {pdf_url}")
                    
                    # Add delay before PDF download
                    time.sleep(random.uniform(1, 3))
                    
                    # Download PDF
                    pdf_response = session.get(pdf_url, timeout=30)
                    if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
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
                            logger.warning(f"PDF text too short ({len(text)} chars)")
                
                # Also try to find direct download links
                for link in soup.find_all('a', href=True):
                    href = link['href']
                    if href.lower().endswith('.pdf'):
                        try:
                            pdf_url = urljoin(domain, href)
                            logger.info(f"Found direct PDF link: {pdf_url}")
                            
                            # Add delay before PDF download
                            time.sleep(random.uniform(1, 3))
                            
                            pdf_response = session.get(pdf_url, timeout=30)
                            if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
                                pdf_file = io.BytesIO(pdf_response.content)
                                pdf_reader = PyPDF2.PdfReader(pdf_file)
                                text = ""
                                for page in pdf_reader.pages:
                                    text += page.extract_text()
                                text = text.strip()
                                if len(text) > 500:
                                    logger.info(f"Got {len(text)} chars from direct PDF")
                                    return text
                                else:
                                    logger.warning(f"Direct PDF text too short ({len(text)} chars)")
                        except Exception as e:
                            logger.error(f"Error with direct PDF link: {e}")
                            continue
                
            except Exception as e:
                logger.error(f"Error with domain {domain}: {e}")
                continue
        
        logger.warning(f"Could not get text from any Sci-Hub domain for {pmid}")
        return None
        
    except Exception as e:
        logger.error(f"Error in Sci-Hub retrieval: {e}")
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

def analyze_literature(drug_name: str) -> pd.DataFrame:
    """Analyze literature for a given drug."""
    try:
        # Get drug names (including synonyms)
        drug_names = get_drug_names(drug_name)
        logger.info(f"Searching for drug names: {drug_names}")
        
        # Search PubMed
        df = search_pubmed_case_reports(drug_name)
        if df.empty:
            logger.warning("No results found in PubMed")
            return pd.DataFrame()
            
        # Get PMIDs
        pmids = df['PMID'].tolist()
        logger.info(f"Found {len(pmids)} papers to analyze")
        
        # Get texts in parallel
        case_reports = get_texts_parallel(pmids)
        logger.info(f"Retrieved {len(case_reports)} full texts")
        
        # Analyze each case report
        for report in case_reports:
            text = report.get('full_text', '')
            if text:
                combined_text = text + ' ' + report['abstract']
                logger.info(f"Text length for {report['pmid']}: {len(combined_text)} characters")
                
                analyzer = CaseReportAnalyzer()
                analyzed = analyzer.analyze(combined_text)
                
                # Update DataFrame with analyzed fields
                idx = df[df['PMID'] == report['pmid']].index[0]
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
                logger.info(f"Found fields for {report['pmid']}: {found_fields}")
        
        # Sort by year descending
        df = df.sort_values('Year', ascending=False)
        
        # Reorder columns
        columns = [
            'Title', 'Abstract', 'PMID', 'DOI', 'PubMed URL', 'Year',
            'Age', 'Sex', 'Oral Dose (mg)', 'QT (ms)', 'QTc (ms)',
            'Heart Rate (bpm)', 'Blood Pressure (mmHg)', 'TdP'
        ]
        df = df[columns]
        
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
