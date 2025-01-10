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
            return float(data["PropertyTable"]["Properties"][0]["MolecularWeight"])
        except Exception as e:
            logger.error(f"Error getting molecular weight: {str(e)}")
            return None

    def _search_chembl_herg(self, compound_name: str) -> Optional[float]:
        """Search ChEMBL for hERG IC50 values"""
        try:
            target = new_client.target
            activity = new_client.activity
            
            # Find hERG target
            herg_targets = target.filter(pref_name__icontains="HERG")
            herg_target_ids = [t['target_chembl_id'] for t in herg_targets]
            
            # Find compound
            molecules = new_client.molecule.filter(pref_name__iexact=compound_name)
            if not molecules:
                logger.warning(f"No molecules found for {compound_name}")
                return None
                
            molecule = molecules[0]
            
            # Get activities
            activities = activity.filter(
                target_chembl_id__in=herg_target_ids,
                molecule_chembl_id=molecule['molecule_chembl_id']
            )
            
            # Find IC50 values
            ic50_values = []
            for act in activities:
                if act.get('standard_type') == 'IC50' and act.get('standard_value'):
                    ic50_values.append(float(act['standard_value']))
            
            if ic50_values:
                # Return median IC50 value
                return sorted(ic50_values)[len(ic50_values)//2]
            
            return None
            
        except Exception as e:
            logger.error(f"Error searching ChEMBL: {str(e)}")
            return None

    def _check_crediblemeds(self, drug_name: str) -> bool:
        """Check if drug is in CredibleMeds list of known TdP risk drugs."""
        try:
            data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
            crediblemeds_file = os.path.join(data_dir, 'crediblemeds_data.txt')
            logger.info(f"Looking for CredibleMeds file at: {crediblemeds_file}")
            logger.info(f"Checking drug: {drug_name}")
            
            if not os.path.exists(crediblemeds_file):
                logger.error(f"CredibleMeds file not found at: {crediblemeds_file}")
                return False
                
            with open(crediblemeds_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if '|' in line:
                        drug, risk = [x.strip() for x in line.split('|')]
                        logger.info(f"Comparing {drug.lower()} with {drug_name.lower()}")
                        if drug.lower() == drug_name.lower():
                            logger.info(f"Found match! Risk: {risk}")
                            return True
            logger.info("No match found in CredibleMeds data")
            return False
        except Exception as e:
            logger.error(f"Error checking CredibleMeds: {str(e)}")
            return False

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
            
            return pd.DataFrame(papers)
        except Exception as e:
            logger.error(f"Error searching literature: {str(e)}")
            return pd.DataFrame()

    def analyze_drug(self, drug_name: str, doses: List[float]) -> Dict:
        """Analyze a drug for hERG interactions"""
        try:
            # Get drug info from PubChem
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            if "IdentifierList" not in data:
                return {"error": f"Drug {drug_name} not found in PubChem"}
            
            cid = data["IdentifierList"]["CID"][0]
            molecular_weight = self._get_molecular_weight(cid)
            
            if not molecular_weight:
                return {"error": "Could not retrieve molecular weight"}
            
            # Get hERG IC50 from ChEMBL
            herg_ic50 = self._search_chembl_herg(drug_name)
            
            # Calculate concentrations
            concentrations = self._calculate_concentrations(
                doses=doses,
                molecular_weight=molecular_weight,
                herg_ic50=herg_ic50
            )
            
            # Check CredibleMeds risk
            is_crediblemeds_risk = self._check_crediblemeds(drug_name)
            
            # Determine if theoretical binding is likely
            theoretical_binding = (
                any(c.ratio_theoretical and c.ratio_theoretical > 0.1 for c in concentrations) 
                if herg_ic50 else None
            ) or is_crediblemeds_risk
            
            # Format concentrations for output
            concentration_data = [
                [c.dose, c.theoretical_max, c.plasma_concentration, c.ratio_theoretical, c.ratio_plasma]
                for c in concentrations
            ]
            
            return {
                "name": drug_name,
                "cid": cid,
                "molecular_weight": molecular_weight,
                "herg_ic50": herg_ic50,
                "herg_source": "ChEMBL Database",
                "concentrations": concentration_data,
                "theoretical_binding": theoretical_binding,
                "crediblemeds_risk": is_crediblemeds_risk,
                "citations": [
                    "Mendez-Lucio O, et al. (2017). Computational Tools for HERG Channel Blockers Safety Assessment in Drug Discovery. Front Pharmacol.",
                    "Gintant G, et al. (2016). Evolution of strategies to improve preclinical cardiac safety testing. Nat Rev Drug Discov.",
                    "Woosley RL, et al. www.Crediblemeds.org, QTdrugs List, AZCERT, Inc. 1822 Innovation Park Dr., Oro Valley, AZ 85755"
                ]
            }
            
        except Exception as e:
            logger.error(f"Error analyzing drug: {str(e)}")
            return {"error": str(e)}
