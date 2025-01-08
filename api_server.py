from flask import Flask, jsonify
from flask_cors import CORS
import os
from dotenv import load_dotenv
import logging
from Bio import Entrez

# Set up logging
logging.basicConfig(level=logging.INFO)
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
    return jsonify({"status": "healthy", "message": "API is running"})

@app.route('/analyze/<drug_name>', methods=['GET'])
def analyze(drug_name):
    try:
        # Check if required environment variables are set
        if not os.environ.get('NCBI_EMAIL') or not os.environ.get('NCBI_API_KEY'):
            logger.error('Missing NCBI credentials')
            return jsonify({'error': 'Missing NCBI credentials'}), 500
            
        logger.info(f'Starting search for drug: {drug_name}')
        
        # Construct search query
        query = f"{drug_name} AND (torsades de pointes OR QT prolongation OR arrhythmia)"
        
        # Search PubMed
        logger.info(f"Executing PubMed search with query: {query}")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            logger.warning(f"No papers found for {drug_name}")
            return jsonify({'message': f"No papers found for {drug_name}"})
            
        # Fetch paper details
        handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
        papers = list(Entrez.parse(handle))
        handle.close()
        
        # Process papers
        results = []
        for paper in papers:
            try:
                paper_info = {
                    'title': paper.get('TI', ''),
                    'authors': '; '.join(paper.get('AU', [])),
                    'year': paper.get('DP', '').split()[0] if paper.get('DP') else None,
                    'journal': paper.get('JT', ''),
                    'abstract': paper.get('AB', ''),
                    'doi': next((id for id in paper.get('AID', []) if id.endswith('[doi]')), '').rstrip('[doi]'),
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
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 10000))
    app.run(host='0.0.0.0', port=port)
