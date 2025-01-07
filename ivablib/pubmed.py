from typing import List, Dict, Union, Optional, Any
from Bio import Entrez, Medline
import requests
from bs4 import BeautifulSoup
import arxiv
from transformers import AutoTokenizer, AutoModel
import torch
import logging
from sentence_transformers import SentenceTransformer
import PyPDF2
from io import BytesIO

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure PubMed API access
Entrez.email = "ghhercock@gmail.com"
Entrez.api_key = "d769d934dc8159bcc9ec9fc29c715a456008"

class FullTextFetcher:
    """Class to handle full text retrieval from various sources"""
    
    def __init__(self):
        """Initialize the fetcher with BioBERT model"""
        try:
            self.tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.2")
            self.model = AutoModel.from_pretrained("dmis-lab/biobert-base-cased-v1.2")
            logger.info("Successfully loaded BioBERT model")
        except Exception as e:
            logger.error(f"Failed to load BioBERT model: {e}")
            raise
    
    def analyze_text(self, text: str) -> Dict[str, any]:
        """
        Analyze text using BioBERT
        
        Args:
            text: Input text to analyze
            
        Returns:
            Dictionary containing analysis results
        """
        try:
            # Truncate text if too long
            max_length = min(self.tokenizer.model_max_length, 512)  # Use smaller chunks
            text = text[:10000]  # Limit initial text size
            
            # Tokenize with truncation
            inputs = self.tokenizer(text, 
                                  return_tensors="pt", 
                                  max_length=max_length, 
                                  truncation=True,
                                  padding=True)
            
            # Get BioBERT embeddings
            with torch.no_grad():
                outputs = self.model(**inputs)
                # Use mean pooling over sequence length
                embeddings = outputs.last_hidden_state.mean(dim=1)
                
                # Convert to Python types for JSON serialization
                embedding_list = embeddings[0].tolist()
                
                return {
                    'text_length': len(text),
                    'embedding_dim': len(embedding_list),
                    'preview': text[:500] + '...' if len(text) > 500 else text,
                    'embedding_sample': embedding_list[:5]  # Just store first 5 values as sample
                }
            
        except Exception as e:
            logger.error(f"Text analysis error: {e}")
            return {
                'error': str(e),
                'text_length': len(text),
                'preview': text[:500] + '...' if len(text) > 500 else text
            }

    def get_from_pmc(self, pmcid: str) -> Optional[str]:
        """Get full text from PubMed Central"""
        try:
            # First try to get PMC ID from PubMed ID
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmcid)
            record = Entrez.read(handle)
            handle.close()
            
            if record[0].get('LinkSetDb') and record[0]['LinkSetDb'][0].get('Link'):
                pmc_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
                logger.info(f"Found PMC ID: {pmc_id}")
                
                # Get full text from PMC
                handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="text", retmode="text")
                full_text = handle.read()
                handle.close()
                
                if full_text:
                    logger.info(f"Retrieved {len(full_text)} chars from PMC")
                    return full_text
            
            # If that fails, try web scraping
            url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmcid}/"
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Try different article content selectors
            article_text = None
            for selector in [
                'div.jig-ncbiinpagenav',
                'article',
                'div.article-body',
                'div.content'
            ]:
                content = soup.find(selector)
                if content:
                    article_text = content.get_text(separator=' ', strip=True)
                    break
            
            if article_text:
                logger.info(f"Retrieved {len(article_text)} chars from PMC web")
                return article_text
            else:
                logger.warning(f"No article text found for PMC{pmcid}")
                return None
                
        except Exception as e:
            logger.error(f"PMC retrieval error: {e}")
            return None

    def get_from_doi(self, doi: str) -> Optional[str]:
        """Get full text using DOI"""
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.5',
                'Connection': 'keep-alive',
            }
            
            # First try direct DOI resolution
            url = f"https://doi.org/{doi}"
            response = requests.get(url, headers=headers, allow_redirects=True)
            response.raise_for_status()
            
            # Check if we got a PDF
            if 'application/pdf' in response.headers.get('Content-Type', '').lower():
                logger.info(f"Retrieved PDF from DOI {doi}")
                return f"PDF available at: {response.url}"
            
            # Try to extract text from HTML
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Try different content selectors
            article_text = None
            for selector in [
                'article',
                'main',
                'div.article-body',
                'div.content',
                'div.fulltext-view',
                'div.article__body',
                'div#article-content'
            ]:
                content = soup.find(selector)
                if content:
                    # Remove unwanted elements
                    for unwanted in content.find_all(['script', 'style', 'nav', 'header', 'footer']):
                        unwanted.decompose()
                    article_text = content.get_text(separator=' ', strip=True)
                    break
            
            if article_text:
                logger.info(f"Retrieved {len(article_text)} chars from DOI {doi}")
                return article_text
            else:
                logger.warning(f"No article text found for DOI {doi}")
                return None
                
        except Exception as e:
            logger.error(f"DOI retrieval error for {doi}: {e}")
            return None

    def get_from_arxiv(self, query: str) -> List[Dict]:
        """Get papers from arXiv"""
        try:
            client = arxiv.Client()
            search = arxiv.Search(query=query, max_results=10)
            
            papers = []
            for result in client.results(search):
                paper_data = {
                    'title': result.title,
                    'abstract': result.summary,
                    'full_text': result.pdf_url,
                    'full_text_source': 'arxiv'
                }
                papers.append(paper_data)
            
            return papers
            
        except Exception as e:
            logger.error(f"arXiv search error: {e}")
            return []

    def get_from_biorxiv(self, query: str) -> List[Dict]:
        """Get papers from bioRxiv/medRxiv"""
        papers = []
        try:
            # Clean up query for bioRxiv API - remove [Title] and other PubMed-specific syntax
            clean_query = query.replace('[Title]', '').strip()
            
            # bioRxiv API requires proper headers
            headers = {
                'User-Agent': 'Mozilla/5.0 (compatible; IvabradineResearch/1.0; +https://github.com/your-repo)',
                'Accept': 'application/json'
            }
            
            # bioRxiv API endpoint - use newer /search endpoint
            # Documentation: https://api.biorxiv.org/
            url = f"https://api.biorxiv.org/search/0/10/{clean_query}"
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            if data.get('collection', []):
                for paper in data['collection']:
                    papers.append({
                        'title': paper.get('title', ''),
                        'doi': paper.get('doi', ''),
                        'year': paper.get('date', '')[:4] if paper.get('date') else '',
                        'abstract': paper.get('abstract', ''),
                        'source': 'biorxiv',
                        'full_text': None
                    })
            
            # medRxiv uses the same API format
            url = f"https://api.medrxiv.org/search/0/10/{clean_query}"
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            if data.get('collection', []):
                for paper in data['collection']:
                    papers.append({
                        'title': paper.get('title', ''),
                        'doi': paper.get('doi', ''),
                        'year': paper.get('date', '')[:4] if paper.get('date') else '',
                        'abstract': paper.get('abstract', ''),
                        'source': 'medrxiv',
                        'full_text': None
                    })
                    
        except Exception as e:
            logger.error(f"bioRxiv/medRxiv retrieval error: {e}")
            
        return papers

    def get_from_scihub(self, doi: str) -> Optional[str]:
        """Get full text from Sci-Hub as a last resort"""
        try:
            # List of Sci-Hub domains (may need updating)
            domains = [
                'https://sci-hub.se',
                'https://sci-hub.st',
                'https://sci-hub.ru'
            ]
            
            for domain in domains:
                try:
                    url = f"{domain}/{doi}"
                    response = requests.get(url, timeout=10)
                    if response.status_code == 200:
                        soup = BeautifulSoup(response.content, 'html.parser')
                        # Look for the PDF iframe or download link
                        pdf_elem = soup.find('iframe', id='pdf') or soup.find('embed', id='pdf')
                        if pdf_elem and pdf_elem.get('src'):
                            pdf_url = pdf_elem['src']
                            if not pdf_url.startswith('http'):
                                pdf_url = domain + pdf_url
                            logger.info(f"Found PDF on Sci-Hub: {pdf_url}")
                            return f"PDF available at: {pdf_url}"
                except Exception as e:
                    continue
                    
        except Exception as e:
            logger.error(f"Sci-Hub retrieval error for {doi}: {e}")
            
        return None

def search_papers(drug_name: str) -> List[Dict[str, Union[str, Dict]]]:
    """
    Search for papers about a drug in PubMed, focusing on case reports and clinical trials.
    Attempts to retrieve full text from PMC, DOI, arXiv, bioRxiv/medRxiv, and Sci-Hub.
    
    Args:
        drug_name: Name of the drug to search for
        
    Returns:
        List of papers with metadata and full text availability
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Searching for papers about {drug_name}...")

    # Load BioBERT model for semantic search
    model = SentenceTransformer('pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb')
    logger.info("Successfully loaded BioBERT model")

    # Construct PubMed query with strict human filter and case reports
    query = f'''((hERG) OR (QT) OR (QTc) OR (torsad*)) AND ({drug_name})) AND 
    ("Humans"[Mesh] NOT (Animals[Mesh] NOT Humans[Mesh])) AND 
    ("Case Reports"[Publication Type] OR "Clinical Trial"[Publication Type])'''
    logger.info(f"Searching PubMed for: {query}")

    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle)
        pmids = record["IdList"]
        logger.info(f"Found {len(pmids)} papers in PubMed")

        papers = []
        for pmid in pmids:
            try:
                # Get paper details from PubMed
                handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                record = Medline.read(handle)
                
                paper = {
                    "title": record.get("TI", "No title"),
                    "source": "pubmed",
                    "pmid": pmid,
                    "doi": record.get("LID", [None])[0] if record.get("LID") else None,
                    "year": record.get("DP", "")[:4],
                    "full_text_available": False,
                    "full_text_source": None
                }

                # Try to get full text from PMC
                try:
                    pmc_id = get_pmc_id(pmid)
                    if pmc_id:
                        text = get_pmc_text(pmc_id)
                        if text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "pmc"
                            continue
                except Exception as e:
                    logger.error(f"PMC retrieval error: {e}")

                # Try DOI if available
                if paper["doi"] and not paper["full_text_available"]:
                    try:
                        text = get_doi_text(paper["doi"])
                        if text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "doi"
                            continue
                    except Exception as e:
                        logger.warning(f"No article text found for DOI {paper['doi']}")

                # Try arXiv search
                if not paper["full_text_available"]:
                    try:
                        arxiv_text = search_arxiv(paper["title"])
                        if arxiv_text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "arxiv"
                            continue
                    except Exception as e:
                        logger.error(f"arXiv search error: {e}")

                # Try bioRxiv/medRxiv search
                if not paper["full_text_available"]:
                    try:
                        # Search bioRxiv
                        biorxiv_text = search_biorxiv(paper["title"])
                        if biorxiv_text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "biorxiv"
                            continue
                        
                        # Search medRxiv
                        medrxiv_text = search_medrxiv(paper["title"])
                        if medrxiv_text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "medrxiv"
                            continue
                    except Exception as e:
                        logger.error(f"bioRxiv/medRxiv search error: {e}")

                # Try Sci-Hub as last resort
                if not paper["full_text_available"] and paper["doi"]:
                    try:
                        scihub_text = get_scihub_text(paper["doi"])
                        if scihub_text:
                            paper["full_text_available"] = True
                            paper["full_text_source"] = "scihub"
                    except Exception as e:
                        pass

                papers.append(paper)
            except Exception as e:
                logger.error(f"Error processing paper {pmid}: {e}")
                continue

        return papers
    except Exception as e:
        logger.error(f"Search error: {e}")
        return []

def search_biorxiv(title: str) -> Optional[str]:
    """Search bioRxiv for a paper by title and return its text if found."""
    try:
        # bioRxiv API endpoint
        url = f"https://api.biorxiv.org/details/biorxiv/{title}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if data["messages"][0]["total"] > 0:
                # Get PDF URL from response
                pdf_url = data["collection"][0]["pdf"]
                if pdf_url:
                    # Download and extract text from PDF
                    response = requests.get(pdf_url)
                    if response.status_code == 200:
                        with BytesIO(response.content) as pdf_file:
                            reader = PyPDF2.PdfReader(pdf_file)
                            text = ""
                            for page in reader.pages:
                                text += page.extract_text()
                            return text
    except Exception as e:
        logger.error(f"bioRxiv search error: {e}")
    return None

def search_medrxiv(title: str) -> Optional[str]:
    """Search medRxiv for a paper by title and return its text if found."""
    try:
        # medRxiv API endpoint
        url = f"https://api.medrxiv.org/details/medrxiv/{title}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if data["messages"][0]["total"] > 0:
                # Get PDF URL from response
                pdf_url = data["collection"][0]["pdf"]
                if pdf_url:
                    # Download and extract text from PDF
                    response = requests.get(pdf_url)
                    if response.status_code == 200:
                        with BytesIO(response.content) as pdf_file:
                            reader = PyPDF2.PdfReader(pdf_file)
                            text = ""
                            for page in reader.pages:
                                text += page.extract_text()
                            return text
    except Exception as e:
        logger.error(f"medRxiv search error: {e}")
    return None

def get_pmc_id(pmid: str) -> Optional[str]:
    """Get PMC ID from PubMed ID"""
    try:
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        record = Entrez.read(handle)
        handle.close()
        
        if record[0].get('LinkSetDb') and record[0]['LinkSetDb'][0].get('Link'):
            pmc_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
            return pmc_id
    except Exception as e:
        logger.error(f"PMC ID retrieval error: {e}")
    return None

def get_pmc_text(pmc_id: str) -> Optional[str]:
    """Get full text from PMC"""
    try:
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="text", retmode="text")
        full_text = handle.read()
        handle.close()
        return full_text
    except Exception as e:
        logger.error(f"PMC retrieval error: {e}")
    return None

def get_doi_text(doi: str) -> Optional[str]:
    """Get full text using DOI"""
    try:
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Connection': 'keep-alive',
        }
        
        # First try direct DOI resolution
        url = f"https://doi.org/{doi}"
        response = requests.get(url, headers=headers, allow_redirects=True)
        response.raise_for_status()
        
        # Check if we got a PDF
        if 'application/pdf' in response.headers.get('Content-Type', '').lower():
            logger.info(f"Retrieved PDF from DOI {doi}")
            return f"PDF available at: {response.url}"
        
        # Try to extract text from HTML
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Try different content selectors
        article_text = None
        for selector in [
            'article',
            'main',
            'div.article-body',
            'div.content',
            'div.fulltext-view',
            'div.article__body',
            'div#article-content'
        ]:
            content = soup.find(selector)
            if content:
                # Remove unwanted elements
                for unwanted in content.find_all(['script', 'style', 'nav', 'header', 'footer']):
                    unwanted.decompose()
                article_text = content.get_text(separator=' ', strip=True)
                break
        
        if article_text:
            logger.info(f"Retrieved {len(article_text)} chars from DOI {doi}")
            return article_text
        else:
            logger.warning(f"No article text found for DOI {doi}")
            return None
            
    except Exception as e:
        logger.error(f"DOI retrieval error for {doi}: {e}")
        return None

def search_arxiv(query: str) -> Optional[str]:
    """Search arXiv for a paper by title and return its text if found."""
    try:
        client = arxiv.Client()
        search = arxiv.Search(query=query, max_results=10)
        
        for result in client.results(search):
            if result.title.lower() == query.lower():
                return result.pdf_url
    except Exception as e:
        logger.error(f"arXiv search error: {e}")
    return None

def get_scihub_text(doi: str) -> Optional[str]:
    """Get full text from Sci-Hub as a last resort"""
    try:
        # List of Sci-Hub domains (may need updating)
        domains = [
            'https://sci-hub.se',
            'https://sci-hub.st',
            'https://sci-hub.ru'
        ]
        
        for domain in domains:
            try:
                url = f"{domain}/{doi}"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    soup = BeautifulSoup(response.content, 'html.parser')
                    # Look for the PDF iframe or download link
                    pdf_elem = soup.find('iframe', id='pdf') or soup.find('embed', id='pdf')
                    if pdf_elem and pdf_elem.get('src'):
                        pdf_url = pdf_elem['src']
                        if not pdf_url.startswith('http'):
                            pdf_url = domain + pdf_url
                        logger.info(f"Found PDF on Sci-Hub: {pdf_url}")
                        return f"PDF available at: {pdf_url}"
            except Exception as e:
                continue
                    
    except Exception as e:
        logger.error(f"Sci-Hub retrieval error for {doi}: {e}")
        
    return None

def main():
    # Test the search
    query = "drug name"
    papers = search_papers(query)
    
    print(f"\nFound {len(papers)} papers")
    for i, paper in enumerate(papers, 1):
        print(f"\nPaper {i}:")
        print(f"Title: {paper.get('title', 'No title')}")
        print(f"PMID: {paper.get('pmid', 'N/A')}")
        print(f"DOI: {paper.get('doi', 'N/A')}")
        
        if paper.get('full_text_available'):
            print(f"Full text retrieved from: {paper.get('full_text_source', 'unknown')}")
        else:
            print("No full text available")
        print("-" * 80)

if __name__ == "__main__":
    main()
