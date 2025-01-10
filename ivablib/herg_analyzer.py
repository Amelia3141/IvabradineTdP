"""hERG Channel Analysis Module"""
import requests
import time
import re
import logging
import os
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime
from chembl_webresource_client.new_client import new_client
import statistics

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

class DrugAnalyzer:
    """Class for analyzing drug properties and hERG interactions"""
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        """Initialize with NCBI credentials"""
        self.email = email
        self.api_key = api_key
    
    def _get_molecular_weight(self, cid: int) -> Optional[float]:
        """Get molecular weight from PubChem"""
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight/JSON"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            logger.info(f"Got molecular weight for CID {cid}: {data['PropertyTable']['Properties'][0]['MolecularWeight']}")
            return float(data["PropertyTable"]["Properties"][0]["MolecularWeight"])
        except Exception as e:
            logger.error(f"Error getting molecular weight: {str(e)}")
            return None

    def _search_chembl_herg(self, compound_name: str) -> Optional[float]:
        """Search ChEMBL for hERG IC50 values"""
        try:
            molecule = new_client.molecule
            activities = new_client.activity
            
            # Search for the molecule
            results = molecule.filter(molecule_synonyms__molecule_synonym__icontains=compound_name).only(['molecule_chembl_id'])
            if not results:
                logger.warning(f"No ChEMBL results found for {compound_name}")
                return None
                
            # Get activities for hERG
            herg_activities = []
            for mol in results:
                acts = activities.filter(
                    molecule_chembl_id=mol['molecule_chembl_id'],
                    target_chembl_id='CHEMBL240',  # hERG
                    standard_type__iexact='IC50'
                ).only(['standard_value', 'standard_units', 'standard_type'])
                
                for act in acts:
                    if 'standard_value' in act and 'standard_units' in act:
                        value = float(act['standard_value'])
                        if act['standard_units'].lower() == 'nm':
                            value = value / 1000  # Convert to μM
                        elif act['standard_units'].lower() == 'um':
                            value = value
                        else:
                            continue
                        herg_activities.append(value)
            
            if not herg_activities:
                logger.warning(f"No hERG activity data found for {compound_name}")
                return None
                
            # Return median IC50
            median_ic50 = statistics.median(herg_activities)
            logger.info(f"Found hERG IC50 for {compound_name}: {median_ic50} μM")
            return median_ic50
            
        except Exception as e:
            logger.error(f"Error searching ChEMBL: {str(e)}")
            return None

    def _check_crediblemeds(self, drug_name: str) -> tuple[bool, str]:
        """Check if drug is in CredibleMeds list of known TdP risk drugs.
        Returns (is_risk, risk_category)"""
        try:
            data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
            crediblemeds_file = os.path.join(data_dir, 'crediblemeds_data.txt')
            logger.info(f"Checking CredibleMeds file: {crediblemeds_file}")
            
            # Create a dictionary of drug names to risk categories for faster lookup
            drug_risks = {}
            with open(crediblemeds_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or '|' not in line:
                        continue
                    drug, risk = [x.strip().lower() for x in line.split('|')]
                    drug_risks[drug] = risk
            
            # Try exact match first
            risk_category = drug_risks.get(drug_name.lower(), '')
            
            # If no exact match, try partial match
            if not risk_category:
                for drug in drug_risks:
                    if drug in drug_name.lower() or drug_name.lower() in drug:
                        risk_category = drug_risks[drug]
                        break
            
            logger.info(f"CredibleMeds result for {drug_name}: {risk_category}")
            return bool(risk_category), risk_category
            
        except Exception as e:
            logger.error(f"Error checking CredibleMeds: {str(e)}")
            return False, ''

    def _calculate_concentrations(self, 
                                doses: List[float], 
                                molecular_weight: float,
                                herg_ic50: Optional[float] = None,
                                volume_of_distribution: float = 100.0,  # L
                                bioavailability: float = 0.4) -> List[DrugConcentration]:
        """Calculate theoretical and plasma concentrations"""
        concentrations = []
        for dose in doses:
            # Convert mg to μmol
            dose_umol = (dose * 1000) / molecular_weight
            
            # Calculate theoretical maximum concentration (μM)
            theoretical_max = dose_umol / volume_of_distribution
            
            # Calculate plasma concentration with bioavailability
            plasma_concentration = theoretical_max * bioavailability
            
            # Calculate ratios if IC50 is available
            ratio_theoretical = theoretical_max / herg_ic50 if herg_ic50 else None
            ratio_plasma = plasma_concentration / herg_ic50 if herg_ic50 else None
            
            concentrations.append(DrugConcentration(
                dose=dose,
                theoretical_max=theoretical_max,
                plasma_concentration=plasma_concentration,
                ratio_theoretical=ratio_theoretical,
                ratio_plasma=ratio_plasma
            ))
        
        logger.info(f"Calculated concentrations: {concentrations}")
        return concentrations

    def _search_literature(self, drug_name: str) -> pd.DataFrame:
        """Search PubChem literature for drug information"""
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/xrefs/PubMedID/JSON"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            papers = []
            if 'InformationList' in data and 'Information' in data['InformationList']:
                pmids = data['InformationList']['Information'][0].get('PubMedID', [])
                
                for pmid in pmids[:10]:  # Limit to first 10 papers
                    try:
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
                        response = requests.get(url)
                        response.raise_for_status()
                        paper_data = response.json()
                        
                        if 'result' in paper_data and str(pmid) in paper_data['result']:
                            paper = paper_data['result'][str(pmid)]
                            papers.append({
                                'Title': paper.get('title', ''),
                                'Authors': ', '.join(paper.get('authors', [])),
                                'Journal': paper.get('source', ''),
                                'Year': paper.get('pubdate', '').split()[0],
                                'PMID': pmid
                            })
                        
                        time.sleep(0.1)  # Rate limiting
                        
                    except Exception as e:
                        logger.warning(f"Error fetching paper {pmid}: {str(e)}")
                        continue
            
            logger.info(f"Found {len(papers)} literature results")
            return pd.DataFrame(papers)
        except Exception as e:
            logger.error(f"Error searching literature: {str(e)}")
            return pd.DataFrame()

    def analyze_drug(self, drug_name: str, doses: List[float]) -> Dict:
        """Analyze a drug for TdP risk"""
        try:
            logger.info(f"Starting analysis for {drug_name}")
            
            # Get hERG data from ChEMBL
            ic50 = self._search_chembl_herg(drug_name)
            logger.info(f"hERG IC50: {ic50}")
            
            # Check CredibleMeds
            is_risk, risk_category = self._check_crediblemeds(drug_name)
            logger.info(f"CredibleMeds: {is_risk}, {risk_category}")
            
            # Get literature
            literature = self._search_literature(drug_name)
            logger.info(f"Found {len(literature)} literature results")
            
            result = {
                'drug_name': drug_name,
                'herg_ic50': ic50,
                'crediblemeds_risk': is_risk,
                'risk_category': risk_category,
                'literature': literature,
                'source': 'ChEMBL Database' if ic50 is not None else 'Unknown'
            }
            
            logger.info(f"Analysis result: {result}")
            return result
            
        except Exception as e:
            logger.error(f"Error analyzing drug: {str(e)}")
            raise
