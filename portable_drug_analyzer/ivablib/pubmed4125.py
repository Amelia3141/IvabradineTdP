import os
import logging
import pandas as pd
import requests
import urllib3
from bs4 import BeautifulSoup
from Bio import Entrez
from Bio import Medline
from typing import List, Dict, Optional
import sys
import time
import re
from urllib.parse import quote
import PyPDF2
import difflib
import csv

# Suppress SSL warnings
urllib3.disable_warnings()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set your email and API key for NCBI from environment variables
logger.info("Setting up NCBI credentials")
Entrez.email = os.environ.get('NCBI_EMAIL', 'your@email.com')
Entrez.api_key = os.environ.get('NCBI_API_KEY')
logger.info(f"Using email: {Entrez.email}")
logger.info("NCBI credentials set up")

def get_from_scihub(url: str, output_dir: str, filename: str) -> Optional[str]:
    """Try to download a paper from Sci-Hub"""
    try:
        # Clean up DOI URL
        clean_url = url.replace('https://doi.org/', '')
        if not clean_url.startswith('http'):
            clean_url = f"https:/doi.org/{clean_url}"
        
        # Try different Sci-Hub domains
        domains = [
            "https://www.sci-hub.se",
            "https://www.sci-hub.ee",
            "https://2024.sci-hub.se",
            "https://dacemirror.sci-hub.se"
        ]
        
        for domain in domains:
            try:
                sci_hub_url = f"{domain}/{clean_url}"
                logger.info(f"Trying URL: {sci_hub_url}")
                
                # Get the Sci-Hub page
                response = requests.get(sci_hub_url, verify=False, timeout=30)
                if response.status_code != 200:
                    continue
                    
                # Parse the page
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Try to find PDF URL in embed tag
                embed = soup.find('embed')
                if embed and 'src' in embed.attrs:
                    pdf_url = embed['src']
                    if pdf_url.startswith('//'):
                        pdf_url = 'https:' + pdf_url
                    elif not pdf_url.startswith('http'):
                        # Add domain to relative URLs
                        if pdf_url.startswith('/'):
                            pdf_url = domain + pdf_url
                        else:
                            # Special handling for dacemirror
                            if 'dacemirror' in pdf_url:
                                pdf_url = 'https://' + pdf_url
                            else:
                                pdf_url = domain + '/' + pdf_url
                    
                    # Remove URL parameters
                    pdf_url = pdf_url.split('#')[0]
                    logger.info(f"Found PDF URL (embed): {pdf_url}")
                    
                    # Download the PDF
                    pdf_path = os.path.join(output_dir, filename)
                    pdf_response = requests.get(pdf_url, verify=False, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
                        return pdf_path
                
                # Try to find PDF URL in iframe
                iframe = soup.find('iframe')
                if iframe and 'src' in iframe.attrs:
                    pdf_url = iframe['src']
                    if pdf_url.startswith('//'):
                        pdf_url = 'https:' + pdf_url
                    elif not pdf_url.startswith('http'):
                        # Add domain to relative URLs
                        if pdf_url.startswith('/'):
                            pdf_url = domain + pdf_url
                        else:
                            # Special handling for dacemirror
                            if 'dacemirror' in pdf_url:
                                pdf_url = 'https://' + pdf_url
                            else:
                                pdf_url = domain + '/' + pdf_url
                    
                    # Remove URL parameters
                    pdf_url = pdf_url.split('#')[0]
                    logger.info(f"Found PDF URL (iframe): {pdf_url}")
                    
                    # Download the PDF
                    pdf_path = os.path.join(output_dir, filename)
                    pdf_response = requests.get(pdf_url, verify=False, timeout=30)
                    if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:  # Ensure it's not an error page
                        with open(pdf_path, 'wb') as f:
                            f.write(pdf_response.content)
                        return pdf_path
                            
            except Exception as e:
                logger.warning(f"Error processing {domain}: {e}")
                continue
                
        return None
        
    except Exception as e:
        logger.error(f"Error downloading from Sci-Hub: {e}")
        return None

def get_doi_from_crossref(title: str) -> Optional[str]:
    """Try to find DOI using CrossRef API with improved title matching"""
    try:
        logger.info(f"Searching CrossRef for: {title}")
        encoded_title = quote(title)
        url = f"https://api.crossref.org/works?query.title={encoded_title}&rows=5&select=DOI,title,score"
        
        headers = {
            'User-Agent': 'IvabradineResearch/1.0 (mailto:your@email.com)',
            'Accept': 'application/json'
        }
        
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()  # Raise exception for non-200 status codes
        
        data = response.json()
        items = data.get('message', {}).get('items', [])
        
        if items:
            # Calculate similarity scores between search title and results
            title_lower = title.lower()
            best_match = None
            best_score = 0
            
            for item in items:
                item_title = item.get('title', [''])[0].lower()
                # Use both CrossRef's score and string similarity
                crossref_score = item.get('score', 0)
                similarity = difflib.SequenceMatcher(None, title_lower, item_title).ratio()
                combined_score = crossref_score * similarity
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_match = item
            
            # Only return DOI if the match is good enough
            if best_match and best_score > 0.5:
                doi = best_match.get('DOI')
                if doi:
                    logger.info(f"Found DOI from CrossRef: {doi} (score: {best_score:.2f})")
                    return doi
            else:
                logger.warning(f"No good match found in CrossRef results (best score: {best_score:.2f})")
                    
    except requests.exceptions.RequestException as e:
        logger.error(f"CrossRef request error for {title}: {e}")
    except ValueError as e:
        logger.error(f"CrossRef JSON parsing error for {title}: {e}")
    except Exception as e:
        logger.error(f"Unexpected CrossRef error for {title}: {e}")
    
    return None

def get_paper_pdf(paper: Dict, output_dir: str) -> Optional[str]:
    """Get PDF for a paper, trying multiple sources"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate safe filename from title
        title = paper.get('title', '')
        safe_title = re.sub(r'[^\w\s-]', '', title)[:100]  # Limit length and remove special chars
        pdf_path = os.path.join(output_dir, f"{safe_title}.pdf")
        
        # Skip if PDF already exists
        if os.path.exists(pdf_path):
            logger.info(f"PDF already exists: {pdf_path}")
            return pdf_path
        
        # Try PMC first if available
        pmc_id = paper.get('pmc')
        if pmc_id:
            logger.info(f"Found PMC ID: {pmc_id}")
            pdf_url = get_pmc_pdf_url(pmc_id)
            if pdf_url:
                if download_pdf(pdf_url, pdf_path):
                    return pdf_path
        
        # Try DOI next
        doi = paper.get('doi')
        if not doi:
            # Try to get DOI from CrossRef
            doi = get_doi_from_crossref(title)
            if doi:
                paper['doi'] = doi
        
        if doi:
            # Try Sci-Hub
            logger.info(f"Trying Sci-Hub for DOI: {doi}")
            pdf_path = get_from_scihub(doi, output_dir, f"{paper.get('pmid', '')}.pdf")
            if pdf_path:
                return pdf_path
        
        logger.warning(f"Could not find PDF for paper: {title}")
        return None
        
    except Exception as e:
        logger.error(f"Error getting PDF for paper {paper.get('title', '')}: {e}")
        return None

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
        '"Humans"[Mesh] AND ("Case Reports"[Publication Type] OR "Cohort Studies"[Publication Type])'
    )
    
    logger.info(f"Searching PubMed: {query}")
    return query

def search_papers(drug_name: str, output_dir: Optional[str] = None) -> List[Dict]:
    """Search for papers about a drug and download PDFs.
    
    Args:
        drug_name: Name of drug to search for
        output_dir: Directory to save PDFs to. If None, will use default.
        
    Returns:
        List of papers with their metadata and PDF paths
    """
    # Set up output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'papers', drug_name.lower())
    os.makedirs(output_dir, exist_ok=True)
    
    # Get drug names from CSV
    drug_names = get_drug_names(drug_name)
    if not drug_names:
        logger.warning(f"No drug names found for {drug_name}")
        return []
        
    logger.info(f"Found drug names: {drug_names}")
    
    # Build and execute query
    query = build_pubmed_query(drug_names)
    logger.info(f"Searching PubMed: {query}")
    
    max_retries = 3
    retry_delay = 1  # seconds
    
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
            record = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            if attempt == max_retries - 1:
                logger.error(f"Failed to search PubMed after {max_retries} attempts: {e}")
                return []
            time.sleep(retry_delay)
            retry_delay *= 2  # exponential backoff
    
    if not record["IdList"]:
        logger.warning("No papers found")
        return []
    
    # Get paper details
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
            papers = list(Medline.parse(handle))
            handle.close()
            break
        except Exception as e:
            if attempt == max_retries - 1:
                logger.error(f"Failed to fetch paper details after {max_retries} attempts: {e}")
                return []
            time.sleep(retry_delay)
            retry_delay *= 2
    
    # Add PDF paths
    results = []
    for paper in papers:
        try:
            # Get DOI
            doi = None
            if 'LID' in paper:
                doi = paper['LID'][0] if paper['LID'][0].endswith('[doi]') else None
            if not doi and 'AID' in paper:
                for aid in paper['AID']:
                    if aid.endswith('[doi]'):
                        doi = aid
                        break
            
            if not doi:
                logger.warning(f"No DOI found for paper: {paper.get('TI', 'Unknown title')}")
                continue
                
            # Clean DOI
            doi = doi.replace('[doi]', '').strip()
            
            # Try to get PDF
            pdf_filename = f"{paper.get('PMID', 'unknown')}.pdf"
            pdf_path = get_from_scihub(doi, output_dir, pdf_filename)
            
            if pdf_path:
                paper['pdf_path'] = pdf_path
                results.append(paper)
            
        except Exception as e:
            logger.error(f"Error processing paper {paper.get('TI', 'Unknown title')}: {e}")
            continue
    
    logger.info(f"Found {len(results)} papers with PDFs")
    return results

def format_pubmed_results(records):
    """Format PubMed results into a DataFrame."""
    results = []
    for article in records:
        try:
            medline = article["MedlineCitation"]
            article_data = medline["Article"]
            
            # Extract authors safely
            authors = article_data.get("AuthorList", [])
            author_names = []
            for author in authors:
                if isinstance(author, dict):
                    last_name = author.get("LastName", "")
                    fore_name = author.get("ForeName", "")
                    author_names.append(f"{last_name} {fore_name}")
            
            # Get publication date
            pub_date = article_data["Journal"]["JournalIssue"]["PubDate"]
            year = pub_date.get("Year", "N/A")
            
            results.append({
                "Title": article_data.get("ArticleTitle", "No title available"),
                "Authors": ", ".join(author_names) if author_names else "No authors listed",
                "Journal": article_data["Journal"].get("Title", "Journal not specified"),
                "Year": year,
                "PMID": medline.get("PMID", "No PMID")
            })
        except Exception as e:
            print(f"Error processing article: {e}")
            continue
            
    return pd.DataFrame(results)

def search_pubmed_case_reports(drug_name):
    """Search PubMed for case reports about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND "case reports"[Publication Type] AND ("torsade de pointes" OR "qt prolongation")'
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Entrez.read(handle)["PubmedArticle"]
        handle.close()

        return format_pubmed_results(records)
    except Exception as e:
        print(f"Error in case reports search: {e}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def search_pubmed_cohort_studies(drug_name):
    """Search PubMed for cohort studies about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND ("cohort studies"[MeSH Terms] OR "cohort"[Title/Abstract]) AND ("torsade de pointes" OR "qt prolongation")'
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Entrez.read(handle)["PubmedArticle"]
        handle.close()

        return format_pubmed_results(records)
    except Exception as e:
        print(f"Error in cohort studies search: {e}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def search_pubmed_clinical_trials(drug_name):
    """Search PubMed for clinical trials about the drug."""
    try:
        query = f'"{drug_name}"[Title/Abstract] AND "clinical trial"[Publication Type] AND ("torsade de pointes" OR "qt prolongation")'
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Entrez.read(handle)["PubmedArticle"]
        handle.close()

        return format_pubmed_results(records)
    except Exception as e:
        print(f"Error in clinical trials search: {e}")
        return pd.DataFrame(columns=["Title", "Authors", "Journal", "Year", "PMID"])

def main():
    if len(sys.argv) != 2:
        print("Usage: python pubmed4125.py <drug_name>")
        sys.exit(1)
        
    drug_name = sys.argv[1].upper()
    papers = search_papers(drug_name)
    print(f"\nFound {len(papers)} papers")

if __name__ == "__main__":
    main()
