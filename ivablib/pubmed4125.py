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
import urllib3

# Suppress SSL warnings
requests.packages.urllib3.disable_warnings()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

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
        logger.info(f"Attempting PMC retrieval for PMID {pmid}")
        pmc_id = get_pmcid_from_pmid(pmid)
        if pmc_id:
            text = get_pmc_text(pmc_id)
            if text:
                source = 'PMC'
                logger.info("Successfully retrieved from PMC")
                return text, source
            else:
                logger.info("No text found in PMC")
    except Exception as e:
        logger.error(f"Error retrieving from PMC: {e}")

    # Try Crossref if we have a DOI
    if doi:
        try:
            logger.info(f"Attempting Crossref retrieval for DOI {doi}")
            text = get_crossref_text(doi)
            if text:
                source = 'Crossref'
                logger.info("Successfully retrieved from Crossref")
                return text, source
            else:
                logger.info("No text found in Crossref or hit paywall")
        except Exception as e:
            logger.error(f"Error retrieving from Crossref: {e}")

    # Finally try Sci-Hub with retries
    max_retries = 3
    for attempt in range(max_retries):
        try:
            logger.info(f"Attempting Sci-Hub retrieval (attempt {attempt + 1}/{max_retries})")
            text = get_scihub_text(pmid, doi)
            if text:
                source = 'Sci-Hub'
                logger.info("Successfully retrieved from Sci-Hub")
                return text, source
            else:
                logger.info(f"No text found in Sci-Hub on attempt {attempt + 1}")
            
            if attempt < max_retries - 1:
                delay = 2 * (attempt + 1)  # Exponential backoff
                logger.info(f"Waiting {delay} seconds before next attempt...")
                time.sleep(delay)
        except Exception as e:
            logger.error(f"Error retrieving from Sci-Hub on attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                delay = 2 * (attempt + 1)
                logger.info(f"Waiting {delay} seconds before next attempt...")
                time.sleep(delay)

    logger.warning(f"Failed to retrieve text from any source for PMID {pmid}")
    return text, source

def get_full_texts(papers):
    """Get full texts for papers"""
    papers_with_text = []
    total_papers = len(papers)
    texts_found = 0
    retrieval_stats = {
        'PMC': 0,
        'Crossref': 0,
        'Sci-Hub': 0,
        'Failed': 0
    }
    
    for paper in papers:
        pmid = str(paper.get('PMID', ''))
        doi = paper.get('DOI', '')
        
        logger.info(f"\nAttempting to retrieve text for paper {pmid}")
        logger.info(f"Title: {paper.get('Title', 'No title')}")
        if doi:
            logger.info(f"DOI: {doi}")
        
        # Get text from various sources
        text, source = get_full_text(pmid, doi)
        
        if text:
            texts_found += 1
            retrieval_stats[source] += 1
            paper_with_text = paper.copy()  # Keep all original paper data
            paper_with_text.update({
                'FullText': text,
                'TextSource': source
            })
            papers_with_text.append(paper_with_text)
            logger.info(f"Successfully retrieved text from {source}")
        else:
            retrieval_stats['Failed'] += 1
            logger.warning(f"Failed to retrieve text for paper {pmid}")
            logger.warning("Attempted sources: PMC, Crossref, and Sci-Hub")
            
    # Log detailed statistics
    logger.info("\nText Retrieval Statistics:")
    logger.info(f"Total papers processed: {total_papers}")
    logger.info(f"Successful retrievals: {texts_found}")
    logger.info(f"Success rate: {(texts_found/total_papers*100):.1f}%")
    logger.info("\nSource breakdown:")
    logger.info(f"- PMC: {retrieval_stats['PMC']} papers")
    logger.info(f"- Crossref: {retrieval_stats['Crossref']} papers")
    logger.info(f"- Sci-Hub: {retrieval_stats['Sci-Hub']} papers")
    logger.info(f"- Failed: {retrieval_stats['Failed']} papers")
    
    return papers_with_text

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
        # First get the redirect URL
        headers = {
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept-Language': 'en-US,en;q=0.9',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
            'Cache-Control': 'max-age=0',
            'Sec-Ch-Ua': '"Not_A Brand";v="8", "Chromium";v="120", "Google Chrome";v="120"',
            'Sec-Ch-Ua-Mobile': '?0',
            'Sec-Ch-Ua-Platform': '"macOS"',
            'Sec-Fetch-Dest': 'document',
            'Sec-Fetch-Mode': 'navigate',
            'Sec-Fetch-Site': 'none',
            'Sec-Fetch-User': '?1'
        }
        
        session = requests.Session()
        session.headers.update(headers)
        
        # Try to get the article page
        response = session.get(f"https://doi.org/{doi}", allow_redirects=True, verify=False, timeout=30)
        final_url = response.url
        
        logger.info(f"Redirected to: {final_url}")
        
        # Get the content from the final URL
        response = session.get(final_url, verify=False, timeout=30)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Check for common paywall indicators
        paywall_indicators = [
            'div[class*="paywall"]',
            'div[id*="paywall"]',
            'div[class*="access-denied"]',
            'div[class*="access-restricted"]',
            'div[class*="subscription-required"]',
            'div[class*="purchase-options"]',
            'div[class*="login-required"]',
            'div[class*="sign-in"]',
            'div[class*="purchase-access"]',
            'div[class*="subscribe"]',
            'div[class*="membership"]'
        ]
        
        for indicator in paywall_indicators:
            if soup.select(indicator):
                logger.info(f"Paywall detected at {final_url}")
                return None
        
        # Try different possible content containers
        selectors = [
            # Common article containers
            ('div', {'class': ['article-body', 'article-text', 'article-content', 'full-text', 'article__body', 'content-main']}),
            ('div', {'class': ['c-article-body', 'c-article-section', 'main-content', 'content-container']}),
            ('div', {'class': ['content-block', 'article-content', 'NLM__article-body', 'main-article']}),
            ('div', {'id': ['body', 'main-content', 'content-main', 'article-content', 'article-body', 'main']}),
            ('article', {}),
            # Publisher-specific containers
            ('div', {'class': ['article-body', 'article__body', 'fulltext', 'content-container']}),  # Wiley, Nature
            ('div', {'class': ['article-body', 'content-block', 'article-content']}),  # BMJ
            ('div', {'class': ['article-body', 'article-text', 'NLM__article-body']}),  # Taylor & Francis
            ('div', {'class': ['article-body', 'article-text', 'article-content', 'full-text']}),  # Elsevier
            ('div', {'class': ['c-article-body', 'c-article-section']}),  # Springer
            # Generic containers
            ('div', {'class': ['content', 'main', 'article', 'paper', 'manuscript']}),
            ('main', {}),
            ('div', {'role': 'main'}),
            ('div', {'id': ['content', 'main', 'article', 'paper']}),
        ]
        
        # Try to find content in any of the selectors
        for tag, attrs in selectors:
            contents = soup.find_all(tag, attrs)
            for content in contents:
                if content:
                    # Remove unwanted elements
                    for unwanted in content.find_all(['div', 'section'], {'class': [
                        'ref-list', 'supplementary-material', 'copyright', 'author-notes',
                        'references', 'footnotes', 'article-footer', 'metrics', 'figure',
                        'acknowledgments', 'peer-review', 'related-content', 'article-comments',
                        'article-tools', 'social-share', 'citation'
                    ]}):
                        unwanted.decompose()
                        
                    # Get text and clean it
                    text = content.get_text(separator=' ', strip=True)
                    text = ' '.join(text.split())  # Normalize whitespace
                    
                    if len(text) > 500:  # Only return if we got substantial text
                        logger.info(f"Successfully extracted {len(text)} characters of text from {final_url}")
                        return text
        
        logger.warning(f"No article text found at {final_url}")
        
    except Exception as e:
        logger.error(f"Error getting Crossref text: {e}")
    return None

def get_scihub_text(pmid, doi=None):
    """Get full text from Sci-Hub."""
    try:
        # Try multiple Sci-Hub domains
        domains = [
            'https://sci-hub.yt',
            'https://sci-hub.is',
            'https://sci-hub.ren',
            'https://sci-hub.wf',
            'https://sci-hub.et',
            'https://sci-hub.st'
        ]
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.9',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        }
        
        session = requests.Session()
        session.headers.update(headers)
        
        for domain in domains:
            try:
                logger.info(f"Trying Sci-Hub domain: {domain}")
                
                # Try DOI first if available, then PMID
                if doi:
                    url = f"{domain}/{doi}"
                else:
                    url = f"{domain}/pubmed/{pmid}"
                
                logger.info(f"Trying URL: {url}")
                response = session.get(url, timeout=30, verify=False, allow_redirects=True)
                
                # Check if we got redirected to the main page (paper not found)
                if response.url.rstrip('/') == domain.rstrip('/'):
                    logger.info(f"Redirected to main page, paper not found on {domain}")
                    continue
                    
                if 'location not found' in response.text.lower() or 'article not found' in response.text.lower():
                    logger.info(f"Paper not found on {domain}")
                    continue
                    
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Look for embedded PDF
                pdf_iframe = soup.find('iframe', {'id': 'pdf'})
                if pdf_iframe and 'src' in pdf_iframe.attrs:
                    pdf_url = pdf_iframe['src']
                    if not pdf_url.startswith('http'):
                        pdf_url = urljoin(domain, pdf_url) if pdf_url.startswith('/') else f"https:{pdf_url}"
                    
                    logger.info(f"Found PDF iframe at {pdf_url}")
                    
                    # Download PDF with retry and exponential backoff
                    for attempt in range(3):
                        try:
                            delay = 2 ** attempt
                            time.sleep(delay)
                            
                            logger.info(f"Downloading PDF from {pdf_url} (attempt {attempt + 1})")
                            pdf_response = session.get(pdf_url, timeout=30, verify=False, allow_redirects=True)
                            
                            # Check content type and size
                            content_type = pdf_response.headers.get('content-type', '').lower()
                            content_length = int(pdf_response.headers.get('content-length', 0))
                            
                            if content_type == 'application/pdf' and content_length > 1000:
                                pdf_file = io.BytesIO(pdf_response.content)
                                try:
                                    pdf_reader = PyPDF2.PdfReader(pdf_file)
                                    text = ""
                                    for page in pdf_reader.pages:
                                        text += page.extract_text()
                                    text = text.strip()
                                    if len(text) > 500:  # Only return if we got substantial text
                                        logger.info(f"Successfully extracted {len(text)} characters from PDF")
                                        return text
                                    else:
                                        logger.warning(f"PDF text too short ({len(text)} chars), might be corrupted")
                                except Exception as e:
                                    logger.error(f"Error reading PDF: {e}")
                            else:
                                logger.warning(f"Invalid PDF response: type={content_type}, size={content_length}")
                        except Exception as e:
                            logger.error(f"Error downloading PDF from {pdf_url} (attempt {attempt + 1}): {e}")
                            if attempt < 2:  # Don't sleep after last attempt
                                continue
                
                # Also try to find direct download links
                download_links = soup.find_all('a', href=True)
                for link in download_links:
                    href = link['href']
                    if href.lower().endswith('.pdf'):
                        try:
                            pdf_url = urljoin(domain, href)
                            logger.info(f"Found direct PDF link: {pdf_url}")
                            
                            pdf_response = session.get(pdf_url, timeout=30, verify=False)
                            if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
                                pdf_file = io.BytesIO(pdf_response.content)
                                pdf_reader = PyPDF2.PdfReader(pdf_file)
                                text = ""
                                for page in pdf_reader.pages:
                                    text += page.extract_text()
                                text = text.strip()
                                if len(text) > 500:
                                    logger.info(f"Successfully extracted {len(text)} characters from PDF")
                                    return text
                                else:
                                    logger.warning(f"PDF text too short ({len(text)} chars), might be corrupted")
                        except Exception as e:
                            logger.error(f"Error downloading PDF from link {pdf_url}: {e}")
                            continue
                
            except Exception as e:
                logger.error(f"Error accessing {domain}: {e}")
                continue
        
        logger.warning(f"Could not retrieve text from any Sci-Hub domain for PMID {pmid}")
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
                    article = paper['MedlineCitation']['Article']
                    
                    # Extract DOI - check multiple possible locations
                    doi = None
                    if 'ELocationID' in article:
                        for eloc in article['ELocationID']:
                            if eloc.attributes.get('EIdType', '') == 'doi':
                                doi = str(eloc)
                                break
                    if not doi and 'ArticleIdList' in paper['PubmedData']:
                        for article_id in paper['PubmedData']['ArticleIdList']:
                            if article_id.attributes.get('IdType', '') == 'doi':
                                doi = str(article_id)
                                break
                    
                    paper_dict = {
                        'PMID': paper['MedlineCitation']['PMID'],
                        'Title': article['ArticleTitle'],
                        'Abstract': article.get('Abstract', {}).get('AbstractText', [''])[0],
                        'Year': article['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                        'DOI': doi
                    }
                    
                    if doi:
                        logger.info(f"Found DOI for paper {paper_dict['PMID']}: {doi}")
                    else:
                        logger.warning(f"No DOI found for paper {paper_dict['PMID']}")
                        
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
            # Get DOI from either ELocationID or ArticleIdList
            doi = None
            if 'ELocationID' in record['MedlineCitation']['Article']:
                for id in record['MedlineCitation']['Article']['ELocationID']:
                    if id.attributes.get('EIdType') == 'doi':
                        doi = str(id)
                        break
            
            if not doi and 'ArticleIdList' in record['PubmedData']:
                for id in record['PubmedData']['ArticleIdList']:
                    if id.attributes.get('IdType') == 'doi':
                        doi = str(id)
                        break

            paper_dict = {
                'PMID': record['MedlineCitation']['PMID'],
                'Title': record['MedlineCitation']['Article']['ArticleTitle'],
                'Abstract': record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                'Year': record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'DOI': doi,
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
        # Wrapper to add delay between requests
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

def process_papers(papers, drug_name):
    """Process papers using CaseReportAnalyzer"""
    from .case_report_analyzer import CaseReportAnalyzer
    
    try:
        analyzer = CaseReportAnalyzer()
        results = analyzer.analyze_papers(papers, drug_name)
        
        if not results.empty:
            return {
                'case_reports': results.to_dict('records'),
                'message': f"Found {len(results)} case reports"
            }
        else:
            return {
                'case_reports': [],
                'message': "No relevant case reports found"
            }
            
    except Exception as e:
        logger.error(f"Error in process_papers: {str(e)}")
        return {
            'error': f"Error processing papers: {str(e)}"
        }

def analyze_literature(drug_name: str):
    """Analyze literature for a given drug."""
    try:
        logger.info(f"Starting literature analysis for {drug_name}")
        
        # Search for papers
        papers_df = search_pubmed_case_reports(drug_name)
        if papers_df.empty:
            logger.info(f"No papers found for {drug_name}")
            return {
                'error': f"No papers found for {drug_name}"
            }
            
        # Convert DataFrame to list of dictionaries
        papers = papers_df.to_dict('records')
        logger.info(f"Found {len(papers)} papers for {drug_name}")
        
        # Get full texts
        papers_with_text = get_full_texts(papers)
        if not papers_with_text:
            return {
                'error': f"No full texts available for {drug_name}. Found {len(papers)} papers but could not retrieve their full text content. This could be because the papers are not freely accessible or are behind a paywall."
            }
            
        logger.info(f"Retrieved {len(papers_with_text)} full texts for {drug_name}")
        
        # Process papers
        results = process_papers(papers_with_text, drug_name)  # Pass drug_name here
        
        # Add summary statistics
        results.update({
            'stats': {
                'total_papers': len(papers),
                'papers_with_text': len(papers_with_text),
                'success_rate': f"{(len(papers_with_text) / len(papers) * 100):.1f}%"
            }
        })
        
        return results
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {str(e)}")
        return {
            'error': f"Error analyzing literature: {str(e)}"
        }

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers.get('case_reports', []))} papers")

if __name__ == "__main__":
    main()
