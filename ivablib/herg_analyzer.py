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
from Bio import Entrez, Medline

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
        self.drug_name_corrections = {
            'donezepil': 'donepezil',
            'galantamin': 'galantamine',
            'rivastigmin': 'rivastigmine',
            'memantine': 'memantine',
            'ivabradin': 'ivabradine'
        }
        
    def _normalize_drug_name(self, drug_name):
        """Normalize drug name to handle common misspellings."""
        drug_name = drug_name.lower().strip()
        # Check for exact matches in corrections
        if drug_name in self.drug_name_corrections:
            return self.drug_name_corrections[drug_name]
        # Check for partial matches (e.g., if missing 'e' at end)
        for misspelling, correct in self.drug_name_corrections.items():
            if drug_name.startswith(misspelling) or misspelling.startswith(drug_name):
                return correct
        return drug_name
        
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

    def _search_chembl_herg(self, compound_name: str) -> Dict:
        """Search ChEMBL for hERG IC50 values."""
        try:
            molecule = new_client.molecule
            target = new_client.target
            activity = new_client.activity
            
            # Get hERG target
            herg_target = 'CHEMBL240'  # This is the ChEMBL ID for hERG
            
            # Search for the drug
            results = molecule.filter(
                molecule_synonyms__molecule_synonym__icontains=compound_name
            ).only(['molecule_chembl_id', 'pref_name'])
            
            logger.info(f"Starting ChEMBL search for {compound_name}")
            
            if not results:
                logger.info(f"No molecules found in ChEMBL for {compound_name}")
                return {'herg_ic50': None, 'source': None}
            
            # Check each molecule for hERG activity
            for mol in results:
                logger.info(f"Checking molecule: {mol['pref_name']} ({mol['molecule_chembl_id']})")
                
                # Get IC50 values for hERG
                activities = activity.filter(
                    molecule_chembl_id=mol['molecule_chembl_id'],
                    target_chembl_id=herg_target,
                    standard_type__iexact='IC50'
                ).only(['standard_value', 'standard_units', 'standard_type'])
                
                ic50_values = []
                for act in activities:
                    if act['standard_value'] and act['standard_units'] == 'nM':
                        # Convert nM to μM
                        ic50_values.append(act['standard_value'] / 1000)
                    elif act['standard_value']:
                        ic50_values.append(act['standard_value'])
                
                if ic50_values:
                    # Use median IC50 value if multiple exist
                    median_ic50 = statistics.median(ic50_values)
                    return {
                        'herg_ic50': median_ic50,
                        'source': f"ChEMBL: {mol['pref_name']} ({mol['molecule_chembl_id']})"
                    }
            
            return {'herg_ic50': None, 'source': None}
            
        except Exception as e:
            logger.error(f"Error searching ChEMBL: {str(e)}")
            return {'herg_ic50': None, 'source': None}

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

    def _calculate_risk_ratio(self, herg_ic50, therapeutic_concentration=None):
        """Calculate the hERG IC50 to therapeutic concentration ratio."""
        if not isinstance(herg_ic50, (int, float)):
            return None
            
        # If no therapeutic concentration provided, use default estimates
        if therapeutic_concentration is None:
            # Default to 1 μM if no specific concentration available
            therapeutic_concentration = 1.0
            
        if therapeutic_concentration > 0:
            return herg_ic50 / therapeutic_concentration
        return None

    def _search_literature(self, drug_name: str) -> pd.DataFrame:
        """Search PubMed for papers about drug and hERG/QT/arrhythmia"""
        try:
            # Set up Entrez
            Entrez.email = self.email
            if self.api_key:
                Entrez.api_key = self.api_key

            # Search PubMed
            query = f"{drug_name} AND (hERG OR QT OR torsade OR arrhythmia)"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                logger.warning("No papers found")
                return pd.DataFrame()

            # Get paper details
            papers = []
            for pmid in record["IdList"]:
                try:
                    time.sleep(0.5)  # Rate limit to avoid 429 errors
                    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                    record = Medline.read(handle)
                    handle.close()

                    if not isinstance(record, dict):
                        continue

                    paper = {
                        'pmid': pmid,
                        'title': record.get('TI', ''),
                        'authors': ', '.join(record.get('AU', [])),  # Join authors into string
                        'journal': record.get('JT', ''),
                        'year': record.get('DP', '').split()[0] if record.get('DP') else '',
                        'abstract': record.get('AB', '')
                    }
                    papers.append(paper)

                except Exception as e:
                    logger.warning(f"Error fetching paper {pmid}: {str(e)}")
                    continue

            logger.info(f"Found {len(papers)} literature results")
            return pd.DataFrame(papers)

        except Exception as e:
            logger.error(f"Error searching literature: {str(e)}")
            return pd.DataFrame()

    def analyze_drug(self, drug_name: str) -> Dict:
        """Analyze a drug for TdP risk"""
        try:
            # Normalize drug name
            drug_name = self._normalize_drug_name(drug_name)
            logger.info(f"Starting analysis for {drug_name}")
            
            # Get hERG data from ChEMBL
            herg_data = self._search_chembl_herg(drug_name)
            herg_ic50 = herg_data['herg_ic50']
            
            # Calculate risk ratio if hERG IC50 is available
            risk_ratio = None
            if isinstance(herg_ic50, (int, float)):
                risk_ratio = self._calculate_risk_ratio(herg_ic50)
            
            # Check CredibleMeds
            is_risk, risk_category = self._check_crediblemeds(drug_name)
            
            # Get literature
            literature = self._search_literature(drug_name)
            
            # Convert None to "Not found" for display
            display_ic50 = "Not found" if herg_ic50 is None else herg_ic50
            
            return {
                'drug_name': drug_name,
                'herg_ic50': display_ic50,
                'herg_source': herg_data['source'] if herg_data['source'] else 'No data available',
                'risk_ratio': risk_ratio,
                'crediblemeds_risk': is_risk,
                'risk_category': risk_category,
                'literature': literature,
                'source': herg_data['source'] if herg_data['source'] else 'No data available'  # Keep for backwards compatibility
            }
            
        except Exception as e:
            logger.error(f"Error analyzing drug: {str(e)}")
            return {"error": str(e)}
