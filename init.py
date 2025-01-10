import sys
print("Python version:", sys.version)

# Import required packages
try:
    import pandas as pd
    print("pandas version:", pd.__version__)
except ImportError:
    print("pandas not installed")

try:
    from Bio import Entrez
    print("biopython version:", Bio.__version__)
except ImportError:
    print("biopython not installed")

try:
    import requests
    print("requests version:", requests.__version__)
except ImportError:
    print("requests not installed")
