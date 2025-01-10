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
    cid: Optional[int]
    molecular_weight: Optional[float]
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
            'User-Agent': 'DrugAnalyzer/1.0',
            'Accept': 'application/json'
        }

    def _make_request(self, url: str, params: Dict = None) -> Optional[requests.Response]:
        """Make API request with rate limiting and error handling"""
        try:
            if params is None:
                params = {}
            
            # Only add API key for PubMed requests
            if 'eutils.ncbi.nlm.nih.gov' in url:
                params['api_key'] = self.api_key
                params['email'] = self.email
            
            response = requests.get(url, headers=self.headers, params=params)
            response.raise_for_status()
            time.sleep(0.34)  # Rate limiting for NCBI APIs
            return response
        except requests.exceptions.RequestException as e:
            if hasattr(e.response, 'status_code') and e.response.status_code == 404:
                return None
            raise Exception(f"API request failed: {str(e)}")

    def get_drug_cid(self, drug_name: str) -> Optional[int]:
        """Get PubChem CID for drug"""
        try:
            url = f"{self.base_url_pubchem}/compound/name/{drug_name}/cids/JSON"
            response = self._make_request(url)
            if response is None:
                return None
            data = response.json()
            return data['IdentifierList']['CID'][0]
        except Exception:
            return None

    def get_molecular_weight(self, cid: int) -> Optional[float]:
        """Get molecular weight for compound"""
        try:
            url = f"{self.base_url_pubchem}/compound/cid/{cid}/property/MolecularWeight/JSON"
            response = self._make_request(url)
            if response is None:
                return None
            data = response.json()
            return data['PropertyTable']['Properties'][0]['MolecularWeight']
        except Exception:
            return None

    def search_herg_data(self, drug_name: str) -> Tuple[Optional[float], Optional[str]]:
        """Search PubMed for hERG IC50 data"""
        try:
            # Search PubMed for relevant articles
            search_url = f"{self.base_url_pubmed}/esearch.fcgi"
            search_params = {
                'db': 'pubmed',
                'term': f"{drug_name} hERG IC50",
                'retmode': 'json',
                'retmax': '5'
            }
            
            response = self._make_request(search_url, search_params)
            if response is None:
                return None, None
                
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
                if response is None:
                    continue
                    
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
            
        except Exception:
            return None, None

    def calculate_concentrations(self, molecular_weight: float, doses: List[float]) -> List[DrugConcentration]:
        """Calculate drug concentrations for given doses"""
        concentrations = []
        for dose in doses:
            # Convert dose from mg to μmol
            dose_in_mol = float(dose) / float(molecular_weight) * 1000.0 if molecular_weight else 0.0
            
            # Calculate concentrations
            theoretical_max = dose_in_mol / 5.0  # Assuming 5L distribution volume
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
            molecular_weight = self.get_molecular_weight(cid) if cid else None
            
            # Get hERG data
            herg_ic50, herg_source = self.search_herg_data(drug_name)
            
            # Calculate concentrations
            concentrations = self.calculate_concentrations(molecular_weight or 0, doses)
            
            # Calculate ratios if hERG IC50 is available
            if herg_ic50:
                for conc in concentrations:
                    try:
                        conc.ratio_theoretical = float(herg_ic50) / float(conc.theoretical_max) if conc.theoretical_max else None
                        conc.ratio_plasma = float(herg_ic50) / float(conc.plasma_concentration) if conc.plasma_concentration else None
                    except (ValueError, TypeError):
                        conc.ratio_theoretical = None
                        conc.ratio_plasma = None
            
            # Determine theoretical binding
            theoretical_binding = False
            if herg_ic50:
                theoretical_binding = any(
                    c.ratio_plasma and c.ratio_plasma < 1 for c in concentrations
                )
            
            # Generate citations
            citations = []
            if cid:
                citations.append(
                    f"National Center for Biotechnology Information (2024). PubChem Compound Summary for CID {cid}, {drug_name}. "
                    f"Retrieved {datetime.now().strftime('%B %d, %Y')}, from https://pubchem.ncbi.nlm.nih.gov/compound/{drug_name}"
                )
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
