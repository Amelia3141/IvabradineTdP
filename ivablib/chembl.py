"""ChEMBL API interaction for hERG data."""
import pandas as pd
from typing import Optional
import requests

def get_herg_ic50(drug_name: str) -> Optional[float]:
    """Get hERG IC50 from CSV."""
    try:
        df = pd.read_csv('hergic50valueschembl.csv')
        drug_data = df[df['Molecule Name'].str.lower() == drug_name.lower()]
        if not drug_data.empty:
            return float(drug_data.iloc[0]['Standard Value'])
    except Exception as e:
        print(f"Error getting hERG IC50: {e}")
    return None

def fetch_pubchem_data(drug_name: str) -> Optional[dict]:
    """Fetch molecular weight and other properties from PubChem."""
    try:
        # Search for compound
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
        response = requests.get(search_url)
        if not response.ok:
            return None
            
        data = response.json()
        if 'IdentifierList' not in data:
            return None
            
        cid = data['IdentifierList']['CID'][0]
        
        # Get properties
        props_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,XLogP/JSON"
        response = requests.get(props_url)
        if not response.ok:
            return None
            
        props = response.json()['PropertyTable']['Properties'][0]
        return {
            'molecular_weight': props.get('MolecularWeight'),
            'logp': props.get('XLogP')
        }
    except Exception as e:
        print(f"PubChem data fetch error: {e}")
        return None
