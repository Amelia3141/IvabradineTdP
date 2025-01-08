from flask import Flask, jsonify, current_app
from flask_cors import CORS
import os
from dotenv import load_dotenv
import logging
from Bio import Entrez, Medline
import time
from Bio import Entrez

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

# Initialize Flask app
app = Flask(__name__)
CORS(app)

# Lazy loading of heavy components
def get_analyzer():
    if not hasattr(current_app, 'analyzer'):
        logger.info("Initializing CaseReportAnalyzer...")
        from ivablib.case_report_analyzer import CaseReportAnalyzer
        current_app.analyzer = CaseReportAnalyzer()
        logger.info("CaseReportAnalyzer initialized successfully")
    return current_app.analyzer

def safe_pubmed_search(query, max_results=5):
    """Safely search PubMed and handle errors"""
    try:
        # First try searching
        logger.info(f"Searching PubMed for: {query}")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        results = Entrez.read(handle, validate=False)
        handle.close()
        
        if not results.get('IdList'):
            logger.warning("No results found")
            return []
            
        # Then fetch the papers
        logger.info(f"Found {len(results['IdList'])} papers, fetching details")
        handle = Entrez.efetch(
            db="pubmed", 
            id=results['IdList'],
            rettype="medline",
            retmode="text"
        )
        
        # Use Medline parser for text format
        records = list(Medline.parse(handle))
        handle.close()
        
        return records
        
    except Exception as e:
        logger.error(f"Error in PubMed search: {str(e)}")
        raise

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
        import Bio
        from Bio import Entrez
        
        # Check environment configuration
        env_status = {
            "NCBI_EMAIL": bool(os.environ.get('NCBI_EMAIL')),
            "NCBI_API_KEY": bool(os.environ.get('NCBI_API_KEY')),
            "FLASK_ENV": os.environ.get('FLASK_ENV', 'production')
        }
        
        # Check if analyzer is initialized
        analyzer_status = "healthy" if hasattr(current_app, 'analyzer') else "not initialized"
        
        # Basic NCBI connection test (only in non-testing env)
        ncbi_status = "skipped"
        if os.environ.get('FLASK_ENV') != 'testing' and env_status["NCBI_EMAIL"] and env_status["NCBI_API_KEY"]:
            try:
                Entrez.email = os.environ.get('NCBI_EMAIL')
                Entrez.api_key = os.environ.get('NCBI_API_KEY')
                handle = Entrez.einfo()
                result = Entrez.read(handle)
                handle.close()
                ncbi_status = "connected" if result else "error"
            except Exception as e:
                ncbi_status = f"error: {str(e)}"
        
        return jsonify({
            "status": "healthy",
            "timestamp": time.time(),
            "environment": env_status,
            "components": {
                "analyzer": analyzer_status,
                "ncbi_api": ncbi_status
            },
            "versions": {
                "python": ".".join(map(str, sys.version_info[:3])),
                "flask": flask.__version__,
                "biopython": Bio.__version__
            }
        })
    except Exception as e:
        return jsonify({
            "status": "unhealthy",
            "error": str(e),
            "timestamp": time.time()
        }), 500

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
        
        # Lazy load the analyzer only when needed
        analyzer = get_analyzer()
        
        # Construct search query
        query = f"{drug_name} AND (torsades de pointes OR QT prolongation OR arrhythmia OR case report)"
        
        # Search PubMed
        papers = safe_pubmed_search(query)
        
        if not papers:
            logger.warning(f"No papers found for {drug_name}")
            return jsonify({
                'drug_name': drug_name,
                'message': f"No papers found for {drug_name}",
                'papers': [],
                'case_reports': [],
                'statistics': {
                    'total_papers': 0,
                    'case_reports_found': 0,
                    'tdp_cases': 0
                }
            })
        
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
        
        # Analyze papers for case reports
        logger.info("Analyzing papers for case reports")
        df = analyzer.analyze_papers(results, drug_name)
        
        # Convert DataFrame to dictionary format
        case_reports = df.to_dict('records') if not df.empty else []
        
        # Calculate statistics
        stats = {
            'total_papers': len(results),
            'case_reports_found': len(case_reports),
            'tdp_cases': len([c for c in case_reports if c.get('had_tdp') == 'Yes'])
        }
        
        logger.info(f"Successfully analyzed {len(results)} papers")
        return jsonify({
            'drug_name': drug_name,
            'statistics': stats,
            'papers': results,
            'case_reports': case_reports
        })
        
    except Exception as e:
        logger.error(f'Error analyzing drug {drug_name}: {str(e)}')
        return jsonify({
            'error': str(e),
            'drug_name': drug_name,
            'environment': {
                'NCBI_EMAIL': bool(os.environ.get('NCBI_EMAIL')),
                'NCBI_API_KEY': bool(os.environ.get('NCBI_API_KEY'))
            }
        }), 500

if __name__ == '__main__':
    import sys
    port = int(os.environ.get('PORT', 10000))
    app.run(host='0.0.0.0', port=port)
