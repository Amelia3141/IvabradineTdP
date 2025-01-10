"""hERG Channel Analysis Module"""
import requests
import time
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime

@dataclass
class DrugConcentration:
    theoretical_max: float  # μM
    plasma_concentration: float  # μM with 40% bioavailability
    ratio_theoretical: Optional[float] = None
    ratio_plasma: Optional[float] = None

@dataclass
class DrugAnalysis:
    name: str
    cid: int
    molecular_weight: float
    herg_ic50: Optional[float]
    herg_source: Optional[str]
    concentrations: List[DrugConcentration]
    theoretical_binding: bool
    citations: List[str]

class DrugAnalyzer:
    def __init__(self, email: str, api_key: str):
        """Initialize with PubMed/NCBI credentials"""
        self.email = email
        self.api_key = api_key
        self.base_url_pubchem = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.base_url_pubmed = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.headers = {
            'User-Agent': f'DrugAnalyzer/1.0 (mailto:{email})'
        }
        self.params = {
            'api_key': api_key,
            'email': email
        }

    def _make_request(self, url: str, params: Dict = None) -> requests.Response:
        """Make API request with rate limiting and error handling"""
        if params:
            params.update(self.params)
        else:
            params = self.params.copy()
        
        response = requests.get(url, headers=self.headers, params=params)
        response.raise_for_status()
        time.sleep(0.34)  # Rate limiting for NCBI APIs (max 3 requests per second)
        return response

    def get_drug_cid(self, drug_name: str) -> int:
        """Get PubChem CID for drug"""
        url = f"{self.base_url_pubchem}/compound/name/{drug_name}/cids/JSON"
        response = self._make_request(url)
        data = response.json()
        return data['IdentifierList']['CID'][0]

    def get_molecular_weight(self, cid: int) -> float:
        """Get molecular weight for compound"""
        url = f"{self.base_url_pubchem}/compound/cid/{cid}/property/MolecularWeight/JSON"
        response = self._make_request(url)
        data = response.json()
        return data['PropertyTable']['Properties'][0]['MolecularWeight']

    def search_herg_data(self, drug_name: str) -> Tuple[Optional[float], Optional[str]]:
        """Search PubMed for hERG IC50 data"""
        # Search PubMed for relevant articles
        search_url = f"{self.base_url_pubmed}/esearch.fcgi"
        search_params = {
            'db': 'pubmed',
            'term': f"{drug_name} hERG IC50",
            'retmode': 'json',
            'retmax': '5'
        }
        
        response = self._make_request(search_url, search_params)
        data = response.json()
        
        if not data['esearchresult']['idlist']:
            return None, None

        # Get abstracts and look for IC50 values
        for pmid in data['esearchresult']['idlist']:
            fetch_url = f"{self.base_url_pubmed}/efetch.fcgi"
            fetch_params = {
                'db': 'pubmed',
                'id': pmid,
                'retmode': 'xml'
            }
            
            response = self._make_request(fetch_url, fetch_params)
            abstract_text = response.text.lower()
            
            # Look for IC50 values in different formats
            import re
            patterns = [
                r'herg.*?ic50.*?(\d+\.?\d*).*?[μµ]m',
                r'ic50.*?(\d+\.?\d*).*?[μµ]m.*?herg',
                r'ic50\s*[=:]\s*(\d+\.?\d*)\s*[μµ]m'
            ]
            
            for pattern in patterns:
                match = re.search(pattern, abstract_text)
                if match:
                    return float(match.group(1)), f"PMID: {pmid}"
        
        return None, None

    def calculate_concentrations(self, molecular_weight: float, doses: List[float]) -> List[DrugConcentration]:
        """Calculate drug concentrations for given doses"""
        concentrations = []
        for dose in doses:
            # Convert dose from mg to μmol
            dose_in_mol = (dose / molecular_weight) * 1000
            
            # Calculate concentrations
            theoretical_max = dose_in_mol / 5  # Assuming 5L distribution volume
            plasma_conc = theoretical_max * 0.4  # 40% bioavailability
            
            concentrations.append(DrugConcentration(
                theoretical_max=theoretical_max,
                plasma_concentration=plasma_conc
            ))
        return concentrations

    def analyze_drug(self, drug_name: str, doses: List[float]) -> DrugAnalysis:
        """Complete drug analysis"""
        try:
            # Get basic drug data
            cid = self.get_drug_cid(drug_name)
            molecular_weight = self.get_molecular_weight(cid)
            
            # Get hERG data
            herg_ic50, herg_source = self.search_herg_data(drug_name)
            
            # Calculate concentrations
            concentrations = self.calculate_concentrations(molecular_weight, doses)
            
            # Calculate ratios if hERG IC50 is available
            if herg_ic50:
                for conc in concentrations:
                    conc.ratio_theoretical = herg_ic50 / conc.theoretical_max
                    conc.ratio_plasma = herg_ic50 / conc.plasma_concentration
            
            # Determine theoretical binding
            theoretical_binding = False
            if herg_ic50:
                theoretical_binding = any(
                    c.ratio_plasma and c.ratio_plasma < 1 for c in concentrations
                )
            
            # Generate citations
            citations = [
                f"National Center for Biotechnology Information (2024). PubChem Compound Summary for CID {cid}, {drug_name}. "
                f"Retrieved {datetime.now().strftime('%B %d, %Y')}, from https://pubchem.ncbi.nlm.nih.gov/compound/{drug_name}"
            ]
            if herg_source:
                citations.append(f"hERG IC50 data source: {herg_source}")
            
            return DrugAnalysis(
                name=drug_name,
                cid=cid,
                molecular_weight=molecular_weight,
                herg_ic50=herg_ic50,
                herg_source=herg_source,
                concentrations=concentrations,
                theoretical_binding=theoretical_binding,
                citations=citations
            )
            
        except Exception as e:
            raise Exception(f"Error analyzing drug {drug_name}: {str(e)}")
