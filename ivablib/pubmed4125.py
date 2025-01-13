import sys
import os
import logging
import time
import random
from typing import List, Dict, Any, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from Bio import Entrez, Medline
import requests
from bs4 import BeautifulSoup
import PyPDF2
import io
import csv
import json

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ivablib.case_report_analyzer import CaseReportAnalyzer

# Suppress SSL warnings
requests.packages.urllib3.disable_warnings()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set up NCBI credentials
Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
Entrez.api_key = os.environ.get('NCBI_API_KEY')
logger.info(f"Using email: {Entrez.email}")
logger.info("NCBI credentials set up")

def get_full_text(pmid: str, doi: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
    """Get full text for a paper using multiple methods."""
    try:
        # Try PMC first
        pmc_id = get_pmcid_from_pmid(pmid)
        if pmc_id:
            logger.info(f"Got PMC ID for {pmid}: {pmc_id}")
            text = get_pmc_text(pmc_id)
            if text and len(text) > 100:  # Ensure we have substantial text
                logger.info(f"Got PMC text for {pmid}: {len(text)} chars")
                return text, 'PMC'
        
        # Try Crossref next if we have a DOI
        if not doi:
            doi = get_doi_from_crossref(pmid)
        if doi:
            logger.info(f"Got DOI for {pmid}: {doi}")
            text = get_crossref_text(doi)
            if text and len(text) > 100:  # Ensure we have substantial text
                logger.info(f"Got Crossref text for {pmid}: {len(text)} chars")
                return text, 'Crossref'
        
        # Try Sci-Hub as last resort
        logger.info(f"Trying Sci-Hub for {pmid}")
        text = get_scihub_text(pmid, doi)
        if text and len(text) > 100:  # Ensure we have substantial text
            logger.info(f"Got Sci-Hub text for {pmid}: {len(text)} chars")
            return text, 'Sci-Hub'
        
        logger.warning(f"Could not get full text for {pmid} from any source")
        return None, None
        
    except Exception as e:
        logger.error(f"Error getting full text for {pmid}: {e}")
        return None, None

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
        url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        }
        
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code != 200:
            logger.warning(f"Got status code {response.status_code} from PMC")
            return None
            
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Get all text content
        text_parts = []
        
        # Get title
        title = soup.find('h1', {'class': 'content-title'})
        if title:
            text_parts.append(title.get_text(strip=True))
        
        # Get abstract
        abstract = soup.find('div', {'class': 'abstract'})
        if abstract:
            text_parts.append(abstract.get_text(strip=True))
        
        # Get main content
        article = soup.find('div', {'class': 'jig-ncbiinpagenav'})
        if not article:
            article = soup.find('article')
        if not article:
            article = soup.find('div', {'class': 'article'})
            
        if article:
            # Get all paragraphs
            for p in article.find_all(['p', 'h2', 'h3', 'h4']):
                text = p.get_text(strip=True)
                if text:  # Only add non-empty paragraphs
                    text_parts.append(text)
        
        if not text_parts:
            logger.warning("No text content found in PMC article")
            return None
            
        # Join all text parts with newlines
        text = '\n\n'.join(text_parts)
        text = ' '.join(text.split())  # Normalize whitespace
        
        if len(text) > 100:  # Only return if we got substantial text
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

def get_scihub_text(pmid: str, doi: Optional[str] = None) -> Optional[str]:
    """Get full text from Sci-Hub."""
    try:
        wait_time = random.uniform(5, 10)
        logger.info(f"Waiting {wait_time:.1f}s before starting Sci-Hub request")
        time.sleep(wait_time)
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5'
        }
        
        # Try different URL patterns
        if doi:
            url = f"https://sci-hub.se/https://doi.org/{doi}"
        else:
            url = f"https://sci-hub.se/pmid/{pmid}"
            
        logger.info(f"Trying Sci-Hub URL: {url}")
        response = requests.get(url, headers=headers, verify=False, timeout=30)
        
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            pdf_element = soup.find('iframe', id='pdf') or soup.find('embed', id='pdf')
            
            if pdf_element and pdf_element.get('src'):
                pdf_url = pdf_element['src']
                if not pdf_url.startswith('http'):
                    pdf_url = f"https:{pdf_url}" if pdf_url.startswith('//') else f"https://sci-hub.se{pdf_url}"
                    
                logger.info(f"Found PDF at {pdf_url}")
                
                # Add delay before PDF request
                wait_time = random.uniform(3, 7)
                time.sleep(wait_time)
                
                pdf_response = requests.get(pdf_url, headers=headers, verify=False, timeout=30)
                if pdf_response.status_code == 200:
                    try:
                        pdf_file = io.BytesIO(pdf_response.content)
                        pdf_reader = PyPDF2.PdfReader(pdf_file)
                        
                        text_parts = []
                        for page in pdf_reader.pages:
                            text = page.extract_text()
                            if text:
                                text_parts.append(text)
                                
                        if text_parts:
                            text = '\n\n'.join(text_parts)
                            if len(text) > 100:
                                logger.info(f"Successfully extracted {len(text)} chars from Sci-Hub PDF")
                                return text
                    except Exception as e:
                        logger.error(f"Error extracting PDF text: {e}")
                        return None
                else:
                    logger.warning(f"Got status code {pdf_response.status_code} for PDF")
        elif response.status_code == 404:
            logger.warning("Paper not found on Sci-Hub")
        elif response.status_code == 400:
            logger.warning("Bad request - might be rate limited or detected as bot")
        else:
            logger.warning(f"Got status code {response.status_code} from Sci-Hub")
            
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
        '('
        'hERG[Title/Abstract] OR '
        'QT[Title/Abstract] OR QTc[Title/Abstract] OR '
        'torsad*[Title/Abstract] OR '
        '"long QT"[Title/Abstract] OR '
        '"prolonged QT"[Title/Abstract] OR '
        '"QT prolongation"[Title/Abstract] OR '
        '"QTc prolongation"[Title/Abstract] OR '
        '"QT interval"[Title/Abstract] OR '
        '"QTc interval"[Title/Abstract] OR '
        '"ventricular tachycardia"[Title/Abstract] OR '
        '"ventricular arrhythmia"[Title/Abstract] OR '
        '"cardiac arrhythmia"[Title/Abstract] OR '
        '"sudden cardiac"[Title/Abstract]'
        ') AND '
        '"Humans"[Mesh] AND ("Case Reports"[Publication Type])'
    )
    
    logger.info(f"Searching PubMed: {query}")
    return query

def search_pubmed_case_reports(drug_name: str) -> pd.DataFrame:
    """Search PubMed for case reports about a drug."""
    try:
        # Initialize Entrez
        Entrez.email = os.getenv('NCBI_EMAIL')
        Entrez.api_key = os.getenv('NCBI_API_KEY')
        
        # Search PubMed
        search_query = f"{drug_name}[Title/Abstract] AND (case report[Title/Abstract] OR case reports[Title/Abstract])"
        handle = Entrez.esearch(db="pubmed", term=search_query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            logger.warning(f"No results found for {drug_name}")
            return pd.DataFrame()
            
        # Get paper details
        handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        papers = []
        
        for record in records:
            try:
                # Get basic paper info
                paper_dict = {
                    'title': record.get('TI', ''),
                    'abstract': record.get('AB', ''),
                    'pmid': record.get('PMID', ''),
                    'doi': '',  # Will be filled by get_doi_from_crossref
                    'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{record.get('PMID', '')}/",
                    'year': record.get('DP', '').split()[0] if record.get('DP') else '',
                    'age': '',
                    'sex': '',
                    'oral_dose_value': '',
                    'oral_dose_unit': '',
                    'oral_dose_freq': '',
                    'qt_value': '',
                    'qtc_value': '',
                    'heart_rate_value': '',
                    'blood_pressure_value': '',
                    'tdp_present': False,
                    'medical_history': '',
                    'medication_history': '',
                    'treatment_course': '',
                    'source': '',
                    'text_length': ''
                }
                papers.append(paper_dict)
            except Exception as e:
                logger.error(f"Error processing paper: {e}")
                continue
                
        handle.close()
        
        # Convert to DataFrame
        df = pd.DataFrame(papers)
        
        # Get DOIs using Crossref
        for idx, row in df.iterrows():
            try:
                if row['pmid']:
                    doi = get_doi_from_crossref(row['pmid'])
                    if doi:
                        df.at[idx, 'doi'] = doi
            except Exception as e:
                logger.error(f"Error getting DOI for {row['pmid']}: {e}")
                continue
        
        return df
        
    except Exception as e:
        logger.error(f"Error in search_pubmed_case_reports: {e}")
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
                'title': record['MedlineCitation']['Article']['ArticleTitle'],
                'abstract': record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0],
                'pmid': pmid,
                'doi': doi or '',
                'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                'year': record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                'age': '',
                'sex': '',
                'oral_dose_value': '',
                'oral_dose_unit': '',
                'oral_dose_freq': '',
                'qt_value': '',
                'qtc_value': '',
                'heart_rate_value': '',
                'blood_pressure_value': '',
                'tdp_present': False,
                'medical_history': '',
                'medication_history': '',
                'treatment_course': '',
                'source': '',
                'text_length': ''
            }
            results.append(paper_dict)
        except Exception as e:
            logger.error(f"Error processing record: {str(e)}")
            continue
            
    return pd.DataFrame(results)

def get_texts_parallel(pmids: List[str]) -> List[Dict[str, Any]]:
    """Get full texts for multiple PMIDs in parallel using ThreadPoolExecutor."""
    try:
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = []
            for pmid in pmids:
                future = executor.submit(get_full_text, pmid)
                futures.append((pmid, future))
            
            texts = []
            for pmid, future in futures:
                try:
                    text, source = future.result(timeout=300)  # 5 minute timeout per paper
                    if text:
                        logger.info(f"Got {len(text)} chars for {pmid}")
                    texts.append({
                        'pmid': pmid,
                        'full_text': text,
                        'source': source,
                        'abstract': ''  # Initialize abstract field
                    })
                except Exception as e:
                    logger.error(f"Error getting text for {pmid}: {e}")
                    texts.append({
                        'pmid': pmid,
                        'full_text': None,
                        'source': None,
                        'abstract': ''
                    })
            
        logger.info(f"Got {len(texts)} full texts out of {len(pmids)} papers")
        return texts
        
    except Exception as e:
        logger.error(f"Error in get_texts_parallel: {e}")
        return []

def analyze_literature(drug_name: str) -> Dict[str, Any]:
    """Analyze literature for a given drug."""
    try:
        # Search PubMed
        df = search_pubmed_case_reports(drug_name)
        if df.empty:
            logger.warning("No results found in PubMed")
            return {'case_reports': []}
            
        # Initialize all required columns
        required_columns = [
            'title',  
            'abstract', 
            'pmid', 
            'doi', 
            'pubmed_url', 
            'year',
            'age',
            'sex',
            'oral_dose_value',
            'oral_dose_unit',
            'oral_dose_freq',
            'qt_value',  
            'qtc_value',  
            'heart_rate_value',
            'blood_pressure_value',
            'tdp_present',  
            'medical_history',  
            'medication_history',  
            'treatment_course',  
            'source',
            'text_length'
        ]
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = ''
                
        # Rename Title to Case Report Title if it exists
        if 'Title' in df.columns:
            df = df.rename(columns={'Title': 'title'})
        
        # Get PMIDs
        pmids = df['pmid'].tolist()
        logger.info(f"Found {len(pmids)} papers to analyze")
        
        # Get texts in parallel
        case_reports = get_texts_parallel(pmids)
        logger.info(f"Retrieved {len(case_reports)} full texts")
        
        # Analyze each case report
        for report in case_reports:
            try:
                pmid = report['pmid']
                text = report.get('full_text', '')
                source = report.get('source', '')
                
                if text:
                    text_length = len(text)
                    logger.info(f"Processing text for {pmid} from {source} ({text_length} chars)")
                    
                    analyzer = CaseReportAnalyzer()
                    analyzed = analyzer.analyze(text)
                    
                    # Update DataFrame with analyzed fields
                    idx = df.index[df['pmid'] == pmid].tolist()
                    if idx:
                        idx = idx[0]
                        # Store the source and text length
                        df.loc[idx, 'source'] = source
                        df.loc[idx, 'text_length'] = text_length
                        
                        # Store analyzed fields
                        if analyzed.get('age'):
                            df.loc[idx, 'age'] = analyzed['age']
                        if analyzed.get('sex'):
                            df.loc[idx, 'sex'] = analyzed['sex']
                        if analyzed.get('oral_dose_value'):
                            df.loc[idx, 'oral_dose_value'] = analyzed['oral_dose_value']
                        if analyzed.get('qt_value'):
                            df.loc[idx, 'qt_value'] = analyzed['qt_value']
                        if analyzed.get('qtc_value'):
                            df.loc[idx, 'qtc_value'] = analyzed['qtc_value']
                        if analyzed.get('heart_rate_value'):
                            df.loc[idx, 'heart_rate_value'] = analyzed['heart_rate_value']
                        if analyzed.get('blood_pressure_value'):
                            df.loc[idx, 'blood_pressure_value'] = analyzed['blood_pressure_value']
                        if analyzed.get('tdp_present') is not None:
                            df.loc[idx, 'tdp_present'] = analyzed['tdp_present']
                        if analyzed.get('medical_history'):
                            df.loc[idx, 'medical_history'] = analyzed['medical_history']
                        if analyzed.get('medication_history'):
                            df.loc[idx, 'medication_history'] = analyzed['medication_history']
                        if analyzed.get('treatment_course'):
                            df.loc[idx, 'treatment_course'] = analyzed['treatment_course']
                        
                        # Log what was found
                        found_fields = {k: v for k, v in analyzed.items() if v}
                        logger.info(f"Found fields for {pmid}: {found_fields}")
                    else:
                        logger.warning(f"Could not find PMID {pmid} in DataFrame")
                else:
                    logger.warning(f"No text available for {pmid}")
                    
            except Exception as e:
                logger.error(f"Error analyzing paper {pmid}: {e}")
                continue
        
        # Sort by year descending
        if 'year' in df.columns:
            df = df.sort_values('year', ascending=False)
        
        # Log summary of papers with data
        data_columns = ['age', 'sex', 'qt_value', 'qtc_value', 'heart_rate_value', 'tdp_present']
        papers_with_data = df[df[data_columns].notna().any(axis=1)]
        logger.info(f"Found data in {len(papers_with_data)} out of {len(df)} papers")
        
        # Log what was found in each paper
        for _, row in papers_with_data.iterrows():
            found_data = {col: row[col] for col in data_columns if pd.notna(row[col]) and row[col] != ''}
            logger.info(f"Paper {row['pmid']} ({row['source']}, {row['text_length']} chars) found: {found_data}")
        
        # Reorder columns to ensure consistent display
        columns = [col for col in required_columns if col in df.columns]
        df = df[columns]
        
        # Convert DataFrame to list of dictionaries
        case_reports = df.to_dict('records')
        
        return {'case_reports': case_reports}
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {e}")
        return {'error': str(e)}

def get_doi_from_crossref(pmid: str) -> Optional[str]:
    """Get DOI for a PMID using Crossref API."""
    try:
        # First try to get DOI from PubMed
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
        records = Entrez.read(handle)
        if records['PubmedArticle']:
            article = records['PubmedArticle'][0]['MedlineCitation']['Article']
            if 'ELocationID' in article:
                for eloc in article['ELocationID']:
                    if eloc.attributes['EIdType'] == 'doi':
                        return str(eloc)
        
        return None
    except Exception as e:
        logger.error(f"Error getting DOI from Crossref for {pmid}: {e}")
        return None

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers['case_reports'])} papers")

if __name__ == "__main__":
    main()
