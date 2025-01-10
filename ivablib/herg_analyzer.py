"""hERG Channel Analysis Module"""
import requests
import time
import re
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from Bio import Entrez, Medline
from datetime import datetime

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(levelname)s:%(name)s:%(message)s'))
    logger.addHandler(handler)

@dataclass
class DrugConcentration:
    dose: float
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
        Entrez.email = email
        Entrez.api_key = api_key

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

    def get_molecular_weight(self, cid: str) -> Optional[float]:
        """Get molecular weight from PubChem"""
        try:
            if not cid:
                logger.warning("No CID provided for molecular weight lookup")
                return None
                
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight/JSON"
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and 'MolecularWeight' in properties[0]:
                        mol_weight = properties[0]['MolecularWeight']
                        try:
                            mol_weight = float(str(mol_weight).strip())
                            if mol_weight > 0:
                                return mol_weight
                            else:
                                logger.warning("Invalid molecular weight: must be positive")
                                return None
                        except (ValueError, TypeError) as e:
                            logger.warning("Could not convert molecular weight to float: {}".format(str(e)))
                            return None
            
            logger.warning("Could not retrieve molecular weight from PubChem")
            return None
            
        except Exception as e:
            logger.error("Error getting molecular weight: {}".format(str(e)))
            return None

    def search_herg_data(self, drug_name: str) -> Tuple[Optional[float], Optional[str]]:
        """Search for hERG IC50 data in literature."""
        try:
            # Search PubMed for articles about the drug and hERG
            query = f"{drug_name} AND (hERG OR KCNH2 OR Kv11.1) AND (IC50 OR IC_50 OR inhibition)"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            
            if not record['IdList']:
                logger.info(f"No articles found for {drug_name}")
                return None, None
                
            # Get details for each article
            handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
            records = Medline.parse(handle)
            records = list(records)  # Convert iterator to list
            handle.close()
            
            # Look for IC50 values in abstracts
            for article in records:
                try:
                    abstract = article.get('AB', '')  # AB is the Medline field for abstract
                    if not abstract:
                        continue
                        
                    # Look for IC50 values with different unit formats
                    ic50_patterns = [
                        r'IC50\s*(?:=|:|\s+of\s+)?\s*(\d+(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?\s*(?:n|μ|micro)?M',
                        r'IC50\s+value\s*(?:of|was|is)?\s*(\d+(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?\s*(?:n|μ|micro)?M',
                        r'hERG.*?(\d+(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?\s*(?:n|μ|micro)?M.*?IC50',
                        r'IC50.*?(\d+(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?\s*(?:n|μ|micro)?M.*?hERG'
                    ]
                    
                    for pattern in ic50_patterns:
                        match = re.search(pattern, abstract, re.IGNORECASE)
                        if match:
                            ic50_value = float(match.group(1))
                            
                            # Get citation
                            authors = article.get('AU', [])
                            first_author = authors[0] if authors else "Unknown"
                            year = article.get('DP', '').split()[0] if article.get('DP') else ''
                            title = article.get('TI', 'Unknown Title')
                            journal = article.get('JT', 'Unknown Journal')
                            
                            citation = f"{first_author} et al. ({year}). {title}. {journal}"
                            logger.info(f"Found hERG IC50 value: {ic50_value} μM")
                            return ic50_value, citation
                            
                except Exception as e:
                    logger.error(f"Error processing article: {str(e)}")
                    continue
            
            logger.info(f"No hERG IC50 values found in literature for {drug_name}")
            return None, None
            
        except Exception as e:
            logger.error(f"Error searching hERG data: {str(e)}")
            return None, None

    def calculate_concentrations(self, molecular_weight: Optional[float], doses: List[float]) -> List[DrugConcentration]:
        """Calculate drug concentrations for given doses"""
        concentrations = []
        try:
            # Check if molecular weight is valid
            if molecular_weight is None:
                logger.warning("No molecular weight available")
                return [DrugConcentration(
                    dose=dose,
                    theoretical_max=None,
                    plasma_concentration=None,
                    ratio_theoretical=None,
                    ratio_plasma=None
                ) for dose in doses]
            
            # Convert molecular weight to float and validate
            try:
                mol_weight = float(str(molecular_weight).strip())
                if mol_weight <= 0:
                    logger.warning("Invalid molecular weight: must be positive")
                    return [DrugConcentration(
                        dose=dose,
                        theoretical_max=None,
                        plasma_concentration=None,
                        ratio_theoretical=None,
                        ratio_plasma=None
                    ) for dose in doses]
                    
                # Calculate concentrations for each dose
                for dose in doses:
                    try:
                        # Convert dose from mg to μmol
                        dose_in_mol = float(str(dose).strip()) / mol_weight * 1000.0
                        
                        # Calculate concentrations
                        theoretical_max = round(dose_in_mol / 5.0, 3)  # Assuming 5L distribution volume
                        plasma_conc = round(theoretical_max * 0.4, 3)  # 40% bioavailability
                        
                        concentrations.append(DrugConcentration(
                            dose=dose,
                            theoretical_max=theoretical_max,
                            plasma_concentration=plasma_conc,
                            ratio_theoretical=None,
                            ratio_plasma=None
                        ))
                    except (ValueError, TypeError, ZeroDivisionError) as e:
                        logger.error("Error calculating concentration for dose {}: {}".format(dose, str(e)))
                        concentrations.append(DrugConcentration(
                            dose=dose,
                            theoretical_max=None,
                            plasma_concentration=None,
                            ratio_theoretical=None,
                            ratio_plasma=None
                        ))
                        
            except (ValueError, TypeError) as e:
                logger.warning("Invalid molecular weight {}: {}".format(molecular_weight, str(e)))
                return [DrugConcentration(
                    dose=dose,
                    theoretical_max=None,
                    plasma_concentration=None,
                    ratio_theoretical=None,
                    ratio_plasma=None
                ) for dose in doses]
                
        except Exception as e:
            logger.error("Error calculating concentrations: {}".format(str(e)))
            concentrations = [DrugConcentration(
                dose=dose,
                theoretical_max=None,
                plasma_concentration=None,
                ratio_theoretical=None,
                ratio_plasma=None
            ) for dose in doses]
            
        return concentrations

    def analyze_drug(self, drug_name: str, doses: List[float]) -> Optional[DrugAnalysis]:
        """Complete drug analysis"""
        try:
            # Get basic drug data
            cid = self.get_drug_cid(drug_name)
            molecular_weight = self.get_molecular_weight(cid) if cid else None
            
            # Get hERG data
            herg_ic50, herg_source = self.search_herg_data(drug_name)
            
            # Calculate concentrations
            concentrations = self.calculate_concentrations(molecular_weight, doses)
            
            # Calculate ratios if hERG IC50 is available
            if herg_ic50 is not None:
                try:
                    herg_ic50 = float(str(herg_ic50).strip())
                    for conc in concentrations:
                        try:
                            if conc.theoretical_max is not None and conc.theoretical_max > 0:
                                conc.ratio_theoretical = round(herg_ic50 / conc.theoretical_max, 3)
                            if conc.plasma_concentration is not None and conc.plasma_concentration > 0:
                                conc.ratio_plasma = round(herg_ic50 / conc.plasma_concentration, 3)
                        except (ValueError, TypeError, ZeroDivisionError) as e:
                            logger.error("Error calculating ratios: {}".format(str(e)))
                            conc.ratio_theoretical = None
                            conc.ratio_plasma = None
                except (ValueError, TypeError) as e:
                    logger.error("Invalid hERG IC50 value {}: {}".format(herg_ic50, str(e)))
            
            # Determine theoretical binding
            theoretical_binding = False
            if herg_ic50 is not None:
                try:
                    theoretical_binding = any(
                        c.ratio_plasma is not None and c.ratio_plasma < 1 
                        for c in concentrations
                    )
                except Exception as e:
                    logger.error("Error determining theoretical binding: {}".format(str(e)))
            
            # Generate citations
            citations = []
            if cid:
                citations.append(
                    "National Center for Biotechnology Information (2024). PubChem Compound Summary for CID {}, {}. "
                    "Retrieved January 10, 2024, from https://pubchem.ncbi.nlm.nih.gov/compound/{}".format(cid, drug_name, drug_name)
                )
            if herg_source:
                citations.append("hERG IC50 data source: {}".format(herg_source))
            
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
            logger.error("Error analyzing drug: {}".format(str(e)))
            return None
