"""hERG Channel Analysis Module"""
import requests
import time
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime
import logging
from Bio import Entrez, Medline
import re

logger = logging.getLogger(__name__)

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
        """Search for hERG IC50 data in literature."""
        try:
            # Search PubMed for articles about the drug and hERG
            query = f"{drug_name} AND (hERG OR KCNH2 OR Kv11.1) AND (IC50 OR IC_50 OR inhibition)"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            
            if not record['IdList']:
                return None, None
                
            # Get details for each article
            handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            # Look for IC50 values in abstracts
            for article in records['PubmedArticle']:
                try:
                    abstract = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [''])[0]
                    if not abstract:
                        continue
                        
                    # Look for IC50 values
                    ic50_pattern = r'IC50\s*(?:=|:|\s+of\s+)?\s*(\d+(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?\s*(?:n|μ|micro)?M'
                    match = re.search(ic50_pattern, abstract, re.IGNORECASE)
                    
                    if match:
                        ic50_value = float(match.group(1))
                        # Get citation
                        authors = article['MedlineCitation']['Article']['AuthorList']
                        first_author = f"{authors[0].get('LastName', '')} {authors[0].get('ForeName', '')}"
                        year = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', '')
                        title = article['MedlineCitation']['Article']['ArticleTitle']
                        journal = article['MedlineCitation']['Article']['Journal']['Title']
                        
                        citation = f"{first_author} et al. ({year}). {title}. {journal}"
                        return ic50_value, citation
                        
                except Exception as e:
                    logger.error(f"Error processing article: {str(e)}")
                    continue
                    
            return None, None
            
        except Exception as e:
            logger.error(f"Error searching hERG data: {str(e)}")
            return None, None

    def calculate_concentrations(self, molecular_weight: float, doses: List[float]) -> List[DrugConcentration]:
        """Calculate drug concentrations for given doses"""
        concentrations = []
        try:
            for dose in doses:
                # Convert dose from mg to μmol
                dose_in_mol = float(dose) / float(molecular_weight) * 1000.0 if molecular_weight else 0.0
                
                # Calculate concentrations
                theoretical_max = dose_in_mol / 5.0  # Assuming 5L distribution volume
                plasma_conc = theoretical_max * 0.4  # 40% bioavailability
                
                concentrations.append(DrugConcentration(
                    theoretical_max=theoretical_max,
                    plasma_concentration=plasma_conc,
                    ratio_theoretical=None,
                    ratio_plasma=None
                ))
        except (ValueError, TypeError) as e:
            logger.error(f"Error calculating concentrations: {str(e)}")
            # Return default concentrations if calculation fails
            concentrations = [DrugConcentration(
                theoretical_max=0.0,
                plasma_concentration=0.0,
                ratio_theoretical=None,
                ratio_plasma=None
            ) for _ in doses]
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
            concentrations = self.calculate_concentrations(molecular_weight or 0.0, doses)
            
            # Calculate ratios if hERG IC50 is available
            if herg_ic50:
                for conc in concentrations:
                    try:
                        if conc.theoretical_max > 0:
                            conc.ratio_theoretical = float(herg_ic50) / float(conc.theoretical_max)
                        if conc.plasma_concentration > 0:
                            conc.ratio_plasma = float(herg_ic50) / float(conc.plasma_concentration)
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
            logger.error(f"Error analyzing drug: {str(e)}")
            return None
