from flask import Flask, jsonify, request
from flask_cors import CORS
from generate_json import analyze_drug
import os
from dotenv import load_dotenv
import logging
import traceback

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

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
            
        logger.info(f'Starting analysis for drug: {drug_name}')
        analysis = analyze_drug(drug_name)
        logger.info(f'Analysis completed for drug: {drug_name}')
        return jsonify(analysis)
    except Exception as e:
        logger.error(f'Error analyzing drug {drug_name}: {str(e)}')
        logger.error(traceback.format_exc())
        return jsonify({'error': str(e), 'details': traceback.format_exc()}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 10000))
    app.run(host='0.0.0.0', port=port)
