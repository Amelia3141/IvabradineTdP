#!/usr/bin/env bash
# exit on error
set -o errexit

echo "Setting up Python environment..."
python -m pip install --upgrade pip

echo "Installing core dependencies first..."
pip install flask==2.3.3 flask-cors==4.0.0 gunicorn==21.2.0 python-dotenv==1.0.0 biopython==1.81 pandas==2.0.3 numpy==1.24.3 requests==2.31.0

echo "Installing PyTorch CPU..."
pip install torch==2.0.1+cpu -f https://download.pytorch.org/whl/cpu

echo "Installing transformer libraries..."
pip install transformers==4.31.0 sentence-transformers==2.2.2

echo "Installing package..."
pip install -e .

echo "Verifying installation..."
python -c "
import sys
import torch
import transformers
import sentence_transformers
print(f'Python version: {sys.version}')
print(f'PyTorch version: {torch.__version__}')
print(f'Transformers version: {transformers.__version__}')
print(f'Sentence-transformers version: {sentence_transformers.__version__}')
"

echo "Build completed successfully"
