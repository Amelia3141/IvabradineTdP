import os
import logging
import pandas as pd
import requests
from bs4 import BeautifulSoup
import PyPDF2
from typing import List, Dict, Optional, Any
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import time
from io import BytesIO
import urllib3
import random
import re
import csv
import sys
from .case_report_analyzer import CaseReportAnalyzer

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Suppress SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email and API key for NCBI from environment variables
logger.info("Setting up NCBI credentials")
from Bio import Entrez
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
        scihub_text = get_text_from_scihub(pmid, doi)
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

def get_crossref_text(doi):
    """Get full text using DOI from Crossref."""
    try:
        # Add proper headers with email
        headers = {
            'User-Agent': 'TdPAnalyzer/1.0 (mailto:your.email@example.com)',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
        }
        
        # Add delay before request
        time.sleep(2)
        
        # Try resolving DOI with proper error handling
        doi_url = f"https://doi.org/{doi}"
        response = requests.get(doi_url, headers=headers, allow_redirects=True, verify=False, timeout=30)
        
        if response.status_code == 403:
            logger.warning(f"Access forbidden (403) for DOI {doi}. Adding delay and retrying...")
            time.sleep(5)  # Longer delay before retry
            response = requests.get(doi_url, headers=headers, allow_redirects=True, verify=False, timeout=30)
        
        if response.status_code == 200:
            final_url = response.url
            logger.info(f"DOI resolved to: {final_url}")
            
            # Add delay before fetching content
            time.sleep(2)
            
            # Get the actual content
            content_response = requests.get(final_url, headers=headers, verify=False, timeout=30)
            if content_response.status_code == 200:
                soup = BeautifulSoup(content_response.text, 'html.parser')
                
                # Try multiple content selectors based on publisher
                if 'sciencedirect.com' in final_url:
                    article = (soup.find('div', {'id': 'centerInner'}) or 
                             soup.find('div', {'id': 'main-content'}) or
                             soup.find('div', {'class': 'article-content'}))
                elif 'springer.com' in final_url:
                    article = (soup.find('div', {'class': 'c-article-body'}) or 
                             soup.find('div', {'id': 'body'}) or
                             soup.find('main', {'class': 'main-content'}))
                elif 'wiley.com' in final_url:
                    article = (soup.find('div', {'class': 'article-content'}) or
                             soup.find('div', {'class': 'main-content'}))
                else:
                    # Generic fallback
                    article = (soup.find('article') or
                             soup.find('div', {'class': ['article', 'content', 'main-content', 'paper']}))
                
                if article:
                    text = article.get_text(separator=' ', strip=True)
                    if len(text) > 1000:  # Sanity check
                        return text
                        
        logger.warning(f"Failed to get text from Crossref for DOI {doi}")
        
    except Exception as e:
        logger.error(f"Error getting Crossref text: {e}")
    return None

def get_text_from_scihub(pmid: str, doi: str) -> Optional[str]:
    """Get full text from Sci-Hub."""
    if not doi:
        return None
        
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
        'Cache-Control': 'max-age=0'
    }
    
    session = requests.Session()
    
    try:
        # First get the main Sci-Hub page
        url = f"https://sci-hub.se/{doi}"
        response = session.get(url, headers=headers, verify=False)
        
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Try to find PDF URL
            pdf_url = None
            
            # Method 1: Check iframe
            iframe = soup.find('iframe', id='pdf')
            if iframe and iframe.get('src'):
                pdf_url = iframe['src']
            
            # Method 2: Check embed
            if not pdf_url:
                embed = soup.find('embed')
                if embed and embed.get('src'):
                    pdf_url = embed['src']
            
            # Method 3: Look for direct PDF link
            if not pdf_url:
                pdf_link = soup.find('a', href=re.compile(r'.*\.pdf$'))
                if pdf_link:
                    pdf_url = pdf_link['href']
            
            if pdf_url:
                if not pdf_url.startswith('http'):
                    pdf_url = f"https:{pdf_url}" if pdf_url.startswith('//') else f"https://sci-hub.se{pdf_url}"
                
                # Get PDF content with same session
                pdf_response = session.get(pdf_url, headers=headers, verify=False)
                if pdf_response.status_code == 200:
                    # Extract text from PDF
                    pdf_file = BytesIO(pdf_response.content)
                    try:
                        reader = PyPDF2.PdfReader(pdf_file)
                        text = ""
                        for page in reader.pages:
                            text += page.extract_text() + "\n"
                        if len(text.strip()) > 100:  # Basic sanity check
                            return text
                    except Exception as e:
                        logger.warning(f"Error extracting text from PDF for {pmid}: {e}")
                        
        elif response.status_code == 403:
            logger.error(f"Access forbidden (403) for {url}. This might be due to rate limiting or IP blocking.")
        else:
            logger.error(f"Failed to access {url}: Status code {response.status_code}")
            
    except Exception as e:
        logger.error(f"Error accessing Sci-Hub for {pmid}: {e}")
        
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

            # Format results using the format_pubmed_results function
            df = format_pubmed_results(records['PubmedArticle'])
            logger.info(f"Processed {len(df)} papers successfully")
            return df
            
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
            pmid = str(record['MedlineCitation']['PMID'])
            
            # Get title and abstract
            article = record['MedlineCitation']['Article']
            title = article['ArticleTitle']
            if hasattr(title, 'attributes'):  # Handle StringElement
                title = str(title)
            elif isinstance(title, (list, str)):
                title = ' '.join(str(t) for t in title) if isinstance(title, list) else str(title)
            
            # Handle abstract text which can be a string or list
            abstract = ''
            if 'Abstract' in article:
                abstract_text = article['Abstract'].get('AbstractText', [])
                if hasattr(abstract_text, 'attributes'):  # Handle StringElement
                    abstract = str(abstract_text)
                elif isinstance(abstract_text, (list, str)):
                    abstract = ' '.join(str(t) for t in abstract_text) if isinstance(abstract_text, list) else str(abstract_text)
            
            # Get DOI
            doi = None
            if 'ArticleIdList' in record['PubmedData']:
                for id_obj in record['PubmedData']['ArticleIdList']:
                    if hasattr(id_obj, 'attributes'):
                        id_type = str(getattr(id_obj.attributes, 'IdType', ''))
                        if id_type.lower() == 'doi':
                            doi = str(id_obj)
                            break
            
            # Get year from PubDate
            year = ''
            pub_date = article['Journal']['JournalIssue']['PubDate']
            if hasattr(pub_date, 'Year'):  # Handle StringElement
                year = str(pub_date.Year)
            elif isinstance(pub_date, dict):
                if 'Year' in pub_date:
                    year = str(pub_date['Year'])
                elif 'MedlineDate' in pub_date:
                    # Try to extract year from MedlineDate
                    medline_date = str(pub_date['MedlineDate'])
                    year_match = re.search(r'\d{4}', medline_date)
                    if year_match:
                        year = year_match.group(0)
            
            # Create result dict with key fields
            paper_dict = {
                'Case Report Title': title,
                'Age': '',
                'Sex': '',
                'Oral Dose (mg)': '',
                'theoretical max concentration (μM)': '',
                '40% bioavailability': '',
                'Theoretical HERG IC50 / Concentration μM': '',
                '40% Plasma concentration': '',
                'Uncorrected QT (ms)': '',
                'QTc': '',
                'QTR': '',
                'QTF': '',
                'Heart Rate (bpm)': '',
                'Torsades de Pointes?': 'No',
                'Blood Pressure (mmHg)': '',
                'Medical History': '',
                'Medication History': '',
                'Course of Treatment': '',
                # Hidden fields for internal use
                'PMID': pmid,
                'DOI': doi or '',
                'PubMed URL': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                'Year': year,
                'Abstract': abstract
            }
            results.append(paper_dict)
        except Exception as e:
            logger.error(f"Error processing record {record.get('MedlineCitation', {}).get('PMID', 'unknown')}: {str(e)}")
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

def analyze_literature(pmids: List[str], drug_name: str) -> pd.DataFrame:
    """Analyze literature for case reports."""
    try:
        # Initialize analyzer
        analyzer = CaseReportAnalyzer()
        
        # Get full texts
        papers = get_texts_parallel(pmids)
        logger.info(f"Retrieved {len(papers)} full texts")
        
        # Analyze each paper
        results = []
        for pmid, text in papers.items():
            try:
                if text:
                    logger.info(f"Text length for {pmid}: {len(text)} characters")
                result = analyzer.analyze_paper({'PMID': pmid, 'FullText': text})
                if result:
                    results.append(result)
                else:
                    logger.warning(f"No results for paper {pmid}")
            except Exception as e:
                logger.error(f"Error analyzing paper {pmid}: {str(e)}")
                continue
        
        # Create DataFrame
        if results:
            df = pd.DataFrame(results)
            logger.info(f"Final dataframe has {len(df)} rows and columns: {list(df.columns)}")
            return df
        else:
            logger.warning("No results to create DataFrame")
            return pd.DataFrame()
            
    except Exception as e:
        logger.error(f"Error in analyze_literature: {str(e)}")
        return pd.DataFrame()

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature([], drug_name)
    print(f"\nFound {len(papers)} papers")

if __name__ == "__main__":
    main()
