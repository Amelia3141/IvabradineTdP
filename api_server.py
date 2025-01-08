from flask import Flask, jsonify
from flask_cors import CORS
import os
from dotenv import load_dotenv
import logging
from Bio import Entrez, Medline
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Configure NCBI credentials
Entrez.email = os.environ.get('NCBI_EMAIL')
Entrez.api_key = os.environ.get('NCBI_API_KEY')

app = Flask(__name__)
CORS(app)

@app.route('/')
def home():
    """Basic health check endpoint"""
    return jsonify({
        "status": "healthy",
        "message": "API is running",
        "timestamp": time.time(),
        "environment": {
            "NCBI_EMAIL": bool(os.environ.get('NCBI_EMAIL')),
            "NCBI_API_KEY": bool(os.environ.get('NCBI_API_KEY'))
        }
    })

@app.route('/health')
def health():
    """Detailed health check endpoint"""
    try:
        # Check if we can import all required modules
        import flask
        import flask_cors
        import dotenv
        import Bio
        
        # Check environment variables
        env_vars = {
            "NCBI_EMAIL": bool(os.environ.get('NCBI_EMAIL')),
            "NCBI_API_KEY": bool(os.environ.get('NCBI_API_KEY')),
            "PORT": os.environ.get('PORT', 10000)
        }
        
        return jsonify({
            "status": "healthy",
            "timestamp": time.time(),
            "environment": env_vars,
            "dependencies": {
                "flask": flask.__version__,
                "flask_cors": flask_cors.__version__,
                "python_dotenv": dotenv.__version__,
                "biopython": Bio.__version__
            }
        })
    except Exception as e:
        logger.error(f"Health check failed: {str(e)}")
        return jsonify({
            "status": "unhealthy",
            "error": str(e)
        }), 500

def safe_pubmed_search(query, max_results=10):
    """Safely search PubMed and handle errors"""
    try:
        # First try searching
        logger.info(f"Searching PubMed for: {query}")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, retmode="xml")
        results = Entrez.read(handle)
        handle.close()
        
        if not results.get('IdList'):
            logger.warning("No results found")
            return []
            
        # Then fetch the papers
        logger.info(f"Found {len(results['IdList'])} papers, fetching details")
        handle = Entrez.efetch(db="pubmed", id=results['IdList'], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        papers = list(records)
        handle.close()
        
        return papers
        
    except Exception as e:
        logger.error(f"Error in PubMed search: {str(e)}")
        raise

@app.route('/analyze/<drug_name>', methods=['GET'])
def analyze(drug_name):
    try:
        # Check if required environment variables are set
        if not os.environ.get('NCBI_EMAIL') or not os.environ.get('NCBI_API_KEY'):
            logger.error('Missing NCBI credentials')
            return jsonify({
                'error': 'Missing NCBI credentials',
                'environment': {
                    'NCBI_EMAIL': bool(os.environ.get('NCBI_EMAIL')),
                    'NCBI_API_KEY': bool(os.environ.get('NCBI_API_KEY'))
                }
            }), 500
            
        logger.info(f'Starting search for drug: {drug_name}')
        
        # Construct search query
        query = f"{drug_name} AND (torsades de pointes OR QT prolongation OR arrhythmia)"
        
        # Search PubMed
        papers = safe_pubmed_search(query)
        
        if not papers:
            logger.warning(f"No papers found for {drug_name}")
            return jsonify({'message': f"No papers found for {drug_name}"})
        
        # Process papers
        results = []
        for paper in papers:
            try:
                paper_info = {
                    'title': paper.get('TI', ''),
                    'authors': '; '.join(paper.get('AU', [])) if paper.get('AU') else '',
                    'year': paper.get('DP', '').split()[0] if paper.get('DP') else None,
                    'journal': paper.get('JT', ''),
                    'abstract': paper.get('AB', ''),
                    'doi': next((id for id in paper.get('AID', []) if '[doi]' in id), '').replace('[doi]', '') if paper.get('AID') else '',
                    'pmid': paper.get('PMID', '')
                }
                results.append(paper_info)
            except Exception as e:
                logger.error(f"Error processing paper: {str(e)}")
                continue
        
        logger.info(f"Successfully found {len(results)} papers")
        return jsonify({
            'drug_name': drug_name,
            'paper_count': len(results),
            'papers': results
        })
        
    except Exception as e:
        logger.error(f'Error analyzing drug {drug_name}: {str(e)}')
        return jsonify({
            'error': str(e),
            'environment': {
                'NCBI_EMAIL': bool(os.environ.get('NCBI_EMAIL')),
                'NCBI_API_KEY': bool(os.environ.get('NCBI_API_KEY'))
            }
        }), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 10000))
    app.run(host='0.0.0.0', port=port)
