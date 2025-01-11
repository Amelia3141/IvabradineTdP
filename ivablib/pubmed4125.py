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
import asyncio
import aiohttp
from functools import lru_cache
import hashlib

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

async def get_full_text_async(pmid: str, session: aiohttp.ClientSession) -> Optional[str]:
    """Get full text of a paper asynchronously"""
    try:
        # Try PMC first
        async with session.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmid}&rettype=txt&retmode=text") as resp:
            if resp.status == 200:
                text = await resp.text()
                if text and len(text) > 100:
                    logger.info(f"Got {pmid} from PMC")
                    return text

        # Get DOI from PubMed
        async with session.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&rettype=xml&retmode=text") as resp:
            if resp.status != 200:
                return None
                
            xml = await resp.text()
            doi_match = re.search(r'<ELocationID EIdType="doi">(.*?)</ELocationID>', xml)
            if not doi_match:
                return None
                
            doi = doi_match.group(1)
            
            # Try Crossref
            async with session.get(f"https://api.crossref.org/works/{doi}") as resp:
                if resp.status == 200:
                    data = await resp.json()
                    if 'message' in data and 'URL' in data['message']:
                        pdf_url = data['message']['URL']
                        try:
                            async with session.get(pdf_url) as pdf_resp:
                                if pdf_resp.status == 200:
                                    pdf_data = await pdf_resp.read()
                                    with open('temp.pdf', 'wb') as f:
                                        f.write(pdf_data)
                                    text = extract_text_from_pdf('temp.pdf')
                                    os.remove('temp.pdf')
                                    if text:
                                        logger.info(f"Got {pmid} from Crossref")
                                        return text
                        except Exception as e:
                            logger.warning(f"Error downloading PDF from Crossref for {pmid}: {str(e)}")
            
            # Finally try Sci-Hub
            text = await get_from_scihub_async(doi, session)
            if text:
                logger.info(f"Got {pmid} from Sci-Hub")
                return text
                
        return None
        
    except Exception as e:
        logger.error(f"Error getting full text for {pmid}: {str(e)}")
        return None

async def get_from_scihub_async(doi: str, session: aiohttp.ClientSession) -> Optional[str]:
    """Get paper from Sci-Hub asynchronously"""
    try:
        scihub_url = f"https://sci-hub.se/{doi}"
        async with session.get(scihub_url) as resp:
            if resp.status != 200:
                return None
                
            html = await resp.text()
            pdf_match = re.search(r'<embed type="application/pdf" src="(.*?)"', html)
            if not pdf_match:
                return None
                
            pdf_url = pdf_match.group(1)
            if not pdf_url.startswith('http'):
                pdf_url = f"https:{pdf_url}"
                
            async with session.get(pdf_url) as pdf_resp:
                if pdf_resp.status != 200:
                    return None
                    
                pdf_data = await pdf_resp.read()
                with open('temp.pdf', 'wb') as f:
                    f.write(pdf_data)
                    
                text = extract_text_from_pdf('temp.pdf')
                os.remove('temp.pdf')
                return text
                
    except Exception as e:
        logger.error(f"Error getting from Sci-Hub: {str(e)}")
        return None

async def get_texts_parallel(pmids: List[str]) -> Dict[str, str]:
    """Get full texts for multiple PMIDs in parallel"""
    async with aiohttp.ClientSession() as session:
        tasks = [get_full_text_async(pmid, session) for pmid in pmids]
        results = await asyncio.gather(*tasks)
        return {pmid: text for pmid, text in zip(pmids, results) if text}

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
        
        # Only search for case reports
        case_reports = search_pubmed_case_reports(drug_name)
        
        # Get full text papers and analyze them
        papers_to_analyze = []
        if not case_reports.empty:
            # Get PMIDs
            pmids = [p for p in case_reports['PMID'].tolist() if p != 'No PMID']
            
            # Get texts in parallel
            texts = asyncio.run(get_texts_parallel(pmids))
            
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
