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

    # Finally try Sci-Hub with retries
    max_retries = 3
    for attempt in range(max_retries):
        try:
            logger.info(f"Trying Sci-Hub for {pmid} (attempt {attempt + 1}/{max_retries})")
            scihub_text = get_scihub_text(pmid, doi)
            if scihub_text:
                text = scihub_text
                source = 'Sci-Hub'
                logger.info(f"Successfully got text from Sci-Hub for {pmid}")
                return text, source
            else:
                logger.info(f"No text found in Sci-Hub for {pmid} on attempt {attempt + 1}")
            # Add a delay between retries
            if attempt < max_retries - 1:
                time.sleep(2)
        except Exception as e:
            logger.error(f"Error getting Sci-Hub text for {pmid} on attempt {attempt + 1}: {e}")
            if attempt < max_retries - 1:
                time.sleep(2)

    logger.warning(f"Could not get full text for {pmid} from any source (PMC, Crossref, or Sci-Hub)")
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
        # First get the redirect URL
        headers = {
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept-Language': 'en-US,en;q=0.9',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        }
        
        response = requests.get(f"https://doi.org/{doi}", headers=headers, allow_redirects=True, verify=False, timeout=30)
        final_url = response.url
        
        logger.info(f"Redirected to: {final_url}")
        
        # Now get the content from the final URL
        response = requests.get(final_url, headers=headers, verify=False, timeout=30)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Check for common paywall indicators
        paywall_indicators = [
            'div[class*="paywall"]',
            'div[id*="paywall"]',
            'div[class*="access-denied"]',
            'div[class*="access-restricted"]',
            'div[class*="subscription-required"]'
        ]
        
        for indicator in paywall_indicators:
            if soup.select(indicator):
                logger.info(f"Paywall detected at {final_url}")
                return None
        
        # Try different possible content containers
        article_text = None
        selectors = [
            # Elsevier
            ('div', {'class': ['article-body', 'article-text', 'article-content', 'full-text']}),
            # Springer
            ('div', {'class': ['c-article-body', 'c-article-section']}),
            # Wiley
            ('div', {'class': ['article-body', 'article__body', 'fulltext']}),
            # Nature
            ('div', {'class': ['c-article-body', 'article-body', 'article__body']}),
            # Generic
            ('div', {'id': ['body', 'main-content', 'content-main']}),
            ('article', {}),
            ('div', {'class': ['content', 'main']}),
        ]
        
        for tag, attrs in selectors:
            content = soup.find(tag, attrs)
            if content:
                # Remove unwanted elements
                for unwanted in content.find_all(['div', 'section'], {'class': ['ref-list', 'supplementary-material', 'copyright', 'author-notes']}):
                    unwanted.decompose()
                article_text = content
                break
        
        if article_text:
            # Clean the text
            text = article_text.get_text(separator=' ', strip=True)
            # Remove excessive whitespace
            text = ' '.join(text.split())
            if len(text) > 500:  # Only return if we got substantial text
                logger.info(f"Successfully extracted {len(text)} characters of text from {final_url}")
                return text
            else:
                logger.info(f"Text too short ({len(text)} chars) from {final_url}, likely paywall or access issue")
                return None
        else:
            logger.warning(f"No article text found at {final_url}")
            
    except Exception as e:
        logger.error(f"Error getting Crossref text: {e}")
    return None

def get_scihub_text(pmid, doi=None):
    """Get full text from Sci-Hub."""
    try:
        # Try multiple Sci-Hub domains
        domains = [
            'https://sci-hub.se',
            'https://sci-hub.st',
            'https://sci-hub.ru',
            'https://sci-hub.ee',
            'https://sci-hub.wf',
            'https://sci-hub.ren',
            'https://sci-hub.cat'
        ]
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.9',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        }
        
        for domain in domains:
            try:
                # Try DOI first if available, then PMID
                if doi:
                    url = f"{domain}/{doi}"
                else:
                    url = f"{domain}/pubmed/{pmid}"
                
                response = requests.get(url, headers=headers, timeout=30, verify=False)
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Look for embedded PDF
                pdf_iframe = soup.find('iframe', {'id': 'pdf'})
                if pdf_iframe and 'src' in pdf_iframe.attrs:
                    pdf_url = pdf_iframe['src']
                    if not pdf_url.startswith('http'):
                        pdf_url = urljoin(domain, pdf_url) if pdf_url.startswith('/') else f"https:{pdf_url}"
                    
                    # Download PDF with retry
                    for _ in range(3):
                        try:
                            pdf_response = requests.get(pdf_url, headers=headers, timeout=30, verify=False)
                            if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
                                pdf_file = io.BytesIO(pdf_response.content)
                                pdf_reader = PyPDF2.PdfReader(pdf_file)
                                text = ""
                                for page in pdf_reader.pages:
                                    text += page.extract_text()
                                if len(text.strip()) > 500:  # Only return if we got substantial text
                                    return text
                                break
                        except Exception as e:
                            logger.error(f"Error downloading PDF from {pdf_url}: {e}")
                            time.sleep(2)
                            continue
                
                # Also try to find direct download links
                download_links = soup.find_all('a', href=True)
                for link in download_links:
                    href = link['href']
                    if href.lower().endswith('.pdf'):
                        try:
                            pdf_url = urljoin(domain, href)
                            pdf_response = requests.get(pdf_url, headers=headers, timeout=30, verify=False)
                            if pdf_response.headers.get('content-type', '').lower() == 'application/pdf':
                                pdf_file = io.BytesIO(pdf_response.content)
                                pdf_reader = PyPDF2.PdfReader(pdf_file)
                                text = ""
                                for page in pdf_reader.pages:
                                    text += page.extract_text()
                                if len(text.strip()) > 500:
                                    return text
                        except Exception as e:
                            logger.error(f"Error downloading PDF from link {pdf_url}: {e}")
                            continue
                
            except Exception as e:
                logger.error(f"Error accessing {domain}: {e}")
                continue
            
            time.sleep(2)  # Add delay between domain attempts
                
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
            paper_with_text = paper.copy()  # Keep all original paper data
            paper_with_text.update({
                'FullText': text,
                'TextSource': source
            })
            papers_with_text.append(paper_with_text)
        else:
            logger.warning(f"No text found for paper {pmid}")
            
    logger.info(f"Found text for {texts_found} out of {total_papers} papers")
    return papers_with_text

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

def process_papers(papers):
    """Process papers using CaseReportAnalyzer"""
    from .case_report_analyzer import CaseReportAnalyzer
    
    try:
        analyzer = CaseReportAnalyzer()
        results = analyzer.analyze_papers(papers)
        
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

def analyze_literature(drug_name: str) -> Dict:
    """Analyze literature for a given drug."""
    try:
        logger.info(f"Starting literature analysis for {drug_name}")
        
        # Search for papers
        papers = search_pubmed_case_reports(drug_name)
        if papers.empty:
            logger.info(f"No papers found for {drug_name}")
            return {
                'case_reports': [],
                'message': f"No case reports found for {drug_name}"
            }
            
        # Convert DataFrame to list of dicts for processing
        papers = papers.to_dict('records')
        logger.info(f"Found {len(papers)} papers to analyze")
        
        # Get full texts
        papers_with_text = get_full_texts(papers)
        
        # Log detailed information about text retrieval
        if papers_with_text:
            texts_by_source = {}
            for paper in papers_with_text:
                source = paper.get('TextSource', 'Unknown')
                texts_by_source[source] = texts_by_source.get(source, 0) + 1
            
            source_breakdown = ', '.join([f"{source}: {count}" for source, count in texts_by_source.items()])
            logger.info(f"Retrieved {len(papers_with_text)} texts: {source_breakdown}")
        
        if not papers_with_text:
            # Create a more informative message about why texts weren't found
            paper_info = []
            for paper in papers:
                pmid = paper.get('PMID', 'Unknown PMID')
                title = paper.get('Title', 'Unknown Title')
                paper_info.append(f"PMID {pmid}: {title}")
            
            paper_details = "\n".join(paper_info)
            logger.info(f"No full texts found for {drug_name}. Papers that were searched:\n{paper_details}")
            
            return {
                'case_reports': [],
                'message': f"No full texts available for {drug_name}. Found {len(papers)} papers but could not retrieve their full text content. This could be because the papers are not freely accessible or are behind a paywall."
            }
        
        # Process papers
        results = process_papers(papers_with_text)
        
        # Add summary information to results
        results['summary'] = {
            'total_papers_found': len(papers),
            'full_texts_retrieved': len(papers_with_text),
            'case_reports_analyzed': len(results.get('case_reports', [])),
            'sources_used': texts_by_source if papers_with_text else {}
        }
        
        return results
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {str(e)}")
        return {
            'error': f"Error analyzing literature: {str(e)}",
            'case_reports': [],
            'message': f"An error occurred while analyzing literature for {drug_name}"
        }

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers['case_reports'])} papers")

if __name__ == "__main__":
    main()
