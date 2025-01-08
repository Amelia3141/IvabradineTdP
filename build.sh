#!/usr/bin/env bash
# exit on error
set -o errexit

# Create cache directory
mkdir -p /opt/render/project/src/.cache/huggingface

# Install Python dependencies
python -m pip install --upgrade pip
pip install -r requirements.txt

# Install the package in editable mode
pip install -e .
