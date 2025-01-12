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
            'Accept': 'text/html',
            'User-Agent': 'TdPAnalyzer/1.0 (mailto:your@email.com)'
        }
        response = requests.get(f"https://doi.org/{doi}", headers=headers, allow_redirects=True)
        final_url = response.url
        
        logger.info(f"Redirected to: {final_url}")
        
        # Now get the content from the final URL
        response = requests.get(final_url, headers=headers)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Try different possible content containers
        article_text = (
            # Elsevier specific selectors
            soup.find('div', {'class': 'article-body'}) or
            soup.find('div', {'id': 'body'}) or
            soup.find('div', {'class': 'article-text'}) or
            soup.find('div', {'class': 'article-content'}) or
            soup.find('div', {'class': 'main-content'}) or
            # ScienceDirect specific
            soup.find('div', {'id': 'centerInner'}) or
            soup.find('div', {'class': 'article'}) or
            # Generic article containers
            soup.find('div', {'class': 'fulltext-view'}) or
            soup.find('article') or
            # Abstract as fallback
            soup.find('div', {'class': 'abstract'}) or
            soup.find('abstract')
        )
        
        if article_text:
            # Clean the text
            text = article_text.get_text(separator=' ', strip=True)
            # Remove excessive whitespace
            text = ' '.join(text.split())
            logger.info(f"Successfully extracted {len(text)} characters of text from {final_url}")
            return text
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
        
        # Get drug names
        drug_names = get_drug_names(drug_name)
        logger.info(f"Found {len(drug_names)} drug names for {drug_name}: {drug_names}")
        
        # Build search query
        drug_terms = ' OR '.join([f'"{name}"[Title/Abstract]' for name in drug_names])
        query = f'({drug_terms}) AND (hERG[Title/Abstract] OR QT[Title/Abstract] OR QTc[Title/Abstract] OR torsad*[Title/Abstract]) AND "Humans"[Mesh] AND ("Case Reports"[Publication Type])'
        
        logger.info(f"Searching PubMed: {query}")
        
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        logger.info("Got esearch handle")
        
        pmids = record["IdList"]
        logger.info(f"Search returned {len(pmids)} results")
        
        if not pmids:
            return {"papers": []}
            
        # Get paper details
        papers = get_paper_details(pmids)
        logger.info(f"Got {len(papers)} paper details")
        
        # Get full texts and extract information
        papers = process_papers(papers)
        logger.info(f"Processed {len(papers)} papers successfully")
        
        # Get full texts
        papers = get_full_texts(papers)
        logger.info(f"Got {len(papers)} full texts out of {len(pmids)} papers")
        
        return {"papers": papers}
        
    except Exception as e:
        logger.error(f"Error in analyze_literature: {str(e)}")
        return {"error": str(e)}

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

def get_full_texts(papers):
    """Get full texts for papers"""
    try:
        # Convert DataFrame to list of dictionaries if it's a DataFrame
        if isinstance(papers, pd.DataFrame):
            papers = papers.to_dict('records')
            
        pmids = [str(paper['PMID']) for paper in papers]
        texts = get_texts_parallel(pmids)
        
        for paper in papers:
            paper['Full Text'] = texts.get(str(paper['PMID']), '')
            
            # If we don't have full text, use the abstract as fallback
            if not paper['Full Text'] and paper.get('Abstract'):
                paper['Full Text'] = paper['Abstract']
                logger.info(f"Using abstract as fallback for PMID {paper['PMID']}")
                
        return papers
        
    except Exception as e:
        logger.error(f"Error getting full texts: {e}")
        return papers

def process_papers(papers):
    """Process papers to extract relevant information"""
    try:
        # Convert DataFrame to list of dictionaries if it's a DataFrame
        if isinstance(papers, pd.DataFrame):
            papers = papers.to_dict('records')
            
        results = []
        for paper in papers:
            pmid = str(paper['PMID'])
            text = paper.get('Full Text', '')
            abstract = paper.get('Abstract', '')
            
            # Combine full text and abstract for better information extraction
            combined_text = f"{text} {abstract}".strip()
            
            report = {
                'PMID': pmid,
                'Title': paper.get('Title', ''),
                'Age': '',
                'Sex': '',
                'Oral Dose (mg)': '',
                'QTc': '',
                'Heart Rate (bpm)': '',
                'Blood Pressure (mmHg)': '',
                'Torsades de Pointes?': 'No',
                'Medical History': '',
                'Medication History': '',
                'Course of Treatment': '',
                'Year': paper.get('Year', ''),
                'Abstract': paper.get('Abstract', ''),
            }
            
            # Extract fields using regex patterns
            if combined_text:
                patterns = {
                    'Age': r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)',
                    'Sex': r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b',
                    'Oral Dose (mg)': r'(?:oral\s+dose|dose[d\s]*orally|given|administered|took|ingested|consumed|prescribed)\s*(?:with|at|of|a|total)?\s*(\d+\.?\d*)\s*(?:mg|milligrams?|g|grams?|mcg|micrograms?)',
                    'QTc': r'(?:QTc[FBR]?|corrected\s+QT|QT\s*corrected|QT[c]?\s*interval\s*(?:corrected)?|corrected\s*QT\s*interval)[\s:]*(?:interval|duration|measurement|prolongation|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?',
                    'Heart Rate (bpm)': r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?',
                    'Blood Pressure (mmHg)': r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:)?\s*([\d/]+)(?:\s*mm\s*Hg)?',
                    'Medical History': r'(?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)[:\.]\s*([^\.]+?)(?:\.|\n|$)',
                    'Medication History': r'(?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)[:\.]\s*([^\.]+?)(?:\.|\n|$)',
                    'Course of Treatment': r'(?:treatment|therapy|management|intervention|administered|given)[:\.]\s*([^\.]+?)(?:\.|\n|$)'
                }
                
                for field, pattern in patterns.items():
                    matches = re.finditer(pattern, combined_text, re.IGNORECASE)
                    for match in matches:
                        if match.groups():
                            # Use the first non-None group
                            value = next((g for g in match.groups() if g is not None), '')
                            if value:
                                report[field] = value
                                break
                
                # Check for TdP specifically
                if re.search(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', combined_text, re.IGNORECASE):
                    report['Torsades de Pointes?'] = 'Yes'
            
            results.append(report)
        
        return results
        
    except Exception as e:
        logger.error(f"Error processing papers: {e}")
        return []

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = analyze_literature(drug_name)
    print(f"\nFound {len(papers['papers'])} papers")

if __name__ == "__main__":
    main()
