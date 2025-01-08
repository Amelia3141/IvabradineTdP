from flask import Flask, jsonify, request
from flask_cors import CORS
from generate_json import analyze_drug
import os

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/')
def home():
    return jsonify({"status": "healthy", "message": "API is running"})

@app.route('/analyze/<drug_name>', methods=['GET'])
def analyze(drug_name):
    try:
        analysis = analyze_drug(drug_name)
        return jsonify(analysis)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 10000))
    app.run(host='0.0.0.0', port=port)
