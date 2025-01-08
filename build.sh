#!/usr/bin/env bash
# exit on error
set -o errexit

# Install Python dependencies
python -m pip install --upgrade pip
pip install -r requirements.txt
