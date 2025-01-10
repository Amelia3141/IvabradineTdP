"""hERG Channel Analysis Module"""
import requests
import time
import re
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
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

class PubChemSearcher:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
    def search_bioassays(self, cid: int) -> List[Dict]:
        """Search PubChem BioAssay database for ion channel data"""
        try:
            # First get list of all bioassays for this compound
            url = f"{self.base_url}/compound/cid/{cid}/aids/JSON"
            response = requests.get(url)
            response.raise_for_status()
            aids = response.json().get('InformationList', {}).get('Information', [{}])[0].get('AID', [])
            
            results = []
            for aid in aids:
                try:
                    # Get assay details
                    url = f"{self.base_url}/assay/aid/{aid}/JSON"
                    response = requests.get(url)
                    response.raise_for_status()
                    assay_data = response.json()
                    
                    # Check if assay is related to ion channels
                    description = str(assay_data).lower()
                    if any(term in description for term in [
                        'herg', 'kcnh2', 'ion channel', 'iks', 'kcnq1', 
                        'calcium channel', 'sodium channel', 'patch clamp'
                    ]):
                        # Get specific results for this compound
                        url = f"{self.base_url}/assay/aid/{aid}/cid/{cid}/JSON"
                        response = requests.get(url)
                        response.raise_for_status()
                        result_data = response.json()
                        
                        # Extract IC50 or other relevant data
                        if 'Table' in result_data:
                            for row in result_data['Table'].get('Row', []):
                                if any(term in str(row).lower() for term in ['ic50', 'inhibition']):
                                    results.append({
                                        'assay_id': aid,
                                        'data': row,
                                        'description': assay_data.get('description', '')
                                    })
                    
                    time.sleep(0.3)  # Rate limiting
                except Exception as e:
                    logger.warning(f"Error processing assay {aid}: {str(e)}")
                    continue
                    
            return results
        except Exception as e:
            logger.error(f"Error searching bioassays: {str(e)}")
            return []

    def search_bioactivity(self, cid: int) -> List[Dict]:
        """Search PubChem's bioactivity data"""
        try:
            url = f"{self.base_url}/compound/cid/{cid}/assaysummary/JSON"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            results = []
            if 'AssaySummaries' in data:
                for assay in data['AssaySummaries']:
                    if any(term in str(assay).lower() for term in [
                        'herg', 'ion channel', 'cardiotoxicity', 'cardiac', 
                        'action potential', 'qt', 'arrhythmia'
                    ]):
                        results.append({
                            'activity_id': assay.get('SID'),
                            'data': assay
                        })
            
            return results
        except Exception as e:
            logger.error(f"Error searching bioactivity: {str(e)}")
            return []

    def get_compound_data(self, cid: int) -> Dict:
        """Get comprehensive compound data from PubChem"""
        try:
            endpoints = [
                ('property', 'MolecularWeight,XLogP,Charge,ComplexityScore,HBondDonorCount,HBondAcceptorCount'),
                ('record', 'synonyms'),
                ('classification', '')
            ]
            
            results = {}
            for endpoint, params in endpoints:
                url = f"{self.base_url}/compound/cid/{cid}/{endpoint}"
                if params:
                    url += f"/{params}"
                url += "/JSON"
                
                response = requests.get(url)
                response.raise_for_status()
                results[endpoint] = response.json()
                time.sleep(0.3)
            
            return results
        except Exception as e:
            logger.error(f"Error getting compound data: {str(e)}")
            return {}

class DrugAnalyzer:
    def __init__(self, email: str, api_key: str):
        """Initialize with PubMed/NCBI credentials"""
        self.email = email
        self.api_key = api_key
        self.pubchem = PubChemSearcher()
        self.base_url_pubchem = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.headers = {
            'User-Agent': 'DrugAnalyzer/1.0',
            'Accept': 'application/json'
        }
        logger.info(f"Initializing DrugAnalyzer with email: {email}")
        if not email or not api_key:
            logger.error("Missing NCBI credentials")
            raise ValueError("NCBI email and API key are required")

    def get_drug_cid(self, drug_name: str) -> Optional[int]:
        """Get PubChem CID for drug"""
        try:
            url = f"{self.base_url_pubchem}/compound/name/{drug_name}/cids/JSON"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            return data['IdentifierList']['CID'][0]
        except Exception as e:
            logger.error(f"Error getting drug CID: {str(e)}")
            return None

    def get_molecular_weight(self, cid: str) -> Optional[float]:
        """Get molecular weight from PubChem"""
        try:
            if not cid:
                logger.warning("No CID provided for molecular weight lookup")
                return None
                
            url = f"{self.base_url_pubchem}/compound/cid/{cid}/property/MolecularWeight/JSON"
            response = requests.get(url)
            response.raise_for_status()
            
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
                            logger.warning(f"Could not convert molecular weight to float: {str(e)}")
                            return None
            
            logger.warning("Could not retrieve molecular weight from PubChem")
            return None
            
        except Exception as e:
            logger.error(f"Error getting molecular weight: {str(e)}")
            return None

    def search_herg_data(self, drug_name: str) -> Tuple[Optional[float], Optional[str]]:
        """Search for hERG IC50 data in PubChem."""
        try:
            # Get CID first
            cid = self.get_drug_cid(drug_name)
            if not cid:
                logger.info(f"No CID found for {drug_name}")
                return None, None
                
            # Search bioassays
            bioassays = self.pubchem.search_bioassays(cid)
            
            # Look for hERG IC50 values in bioassays
            for assay in bioassays:
                data_str = str(assay['data']).lower()
                if 'herg' in data_str or 'kcnh2' in data_str:
                    # Look for IC50 values
                    ic50_pattern = r'ic50\s*[=:]\s*(\d+(?:\.\d+)?)\s*(?:n|µ|u|micro|m|nano|μ)?m'
                    match = re.search(ic50_pattern, data_str)
                    if match:
                        try:
                            ic50_value = float(match.group(1))
                            if ic50_value > 0:
                                citation = f"PubChem BioAssay AID: {assay['assay_id']}"
                                logger.info(f"Found hERG IC50 value: {ic50_value} μM")
                                return ic50_value, citation
                        except (ValueError, TypeError) as e:
                            logger.warning(f"Error converting IC50 value: {str(e)}")
                            continue
            
            # If no IC50 found in bioassays, check bioactivity data
            bioactivities = self.pubchem.search_bioactivity(cid)
            for activity in bioactivities:
                data_str = str(activity['data']).lower()
                if 'herg' in data_str or 'kcnh2' in data_str:
                    ic50_pattern = r'ic50\s*[=:]\s*(\d+(?:\.\d+)?)\s*(?:n|µ|u|micro|m|nano|μ)?m'
                    match = re.search(ic50_pattern, data_str)
                    if match:
                        try:
                            ic50_value = float(match.group(1))
                            if ic50_value > 0:
                                citation = f"PubChem BioActivity SID: {activity['activity_id']}"
                                logger.info(f"Found hERG IC50 value: {ic50_value} μM")
                                return ic50_value, citation
                        except (ValueError, TypeError) as e:
                            logger.warning(f"Error converting IC50 value: {str(e)}")
                            continue
            
            logger.info(f"No hERG IC50 values found in PubChem for {drug_name}")
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
                        dose_umol = (dose * 1000) / mol_weight
                        
                        # Estimate volume of distribution (L)
                        volume = 70  # Average adult body weight in kg
                        
                        # Calculate theoretical maximum concentration (μM)
                        theoretical_max = dose_umol / volume
                        
                        # Calculate plasma concentration with 40% bioavailability
                        plasma_concentration = theoretical_max * 0.4
                        
                        concentrations.append(DrugConcentration(
                            dose=dose,
                            theoretical_max=theoretical_max,
                            plasma_concentration=plasma_concentration,
                            ratio_theoretical=None,
                            ratio_plasma=None
                        ))
                        
                    except Exception as e:
                        logger.error(f"Error calculating concentration for dose {dose}: {str(e)}")
                        concentrations.append(DrugConcentration(
                            dose=dose,
                            theoretical_max=None,
                            plasma_concentration=None,
                            ratio_theoretical=None,
                            ratio_plasma=None
                        ))
                        
            except (ValueError, TypeError) as e:
                logger.error(f"Error converting molecular weight: {str(e)}")
                return [DrugConcentration(
                    dose=dose,
                    theoretical_max=None,
                    plasma_concentration=None,
                    ratio_theoretical=None,
                    ratio_plasma=None
                ) for dose in doses]
                
            return concentrations
            
        except Exception as e:
            logger.error(f"Error calculating concentrations: {str(e)}")
            return [DrugConcentration(
                dose=dose,
                theoretical_max=None,
                plasma_concentration=None,
                ratio_theoretical=None,
                ratio_plasma=None
            ) for dose in doses]

    def analyze_drug(self, drug_name: str, doses: List[float]) -> Dict:
        """Complete drug analysis"""
        # Initialize default response
        response = {
            'name': drug_name,
            'cid': None,
            'molecular_weight': None,
            'herg_ic50': None,
            'herg_source': None,
            'concentrations': [],
            'theoretical_binding': False,
            'citations': []
        }
        
        try:
            # Get basic drug data
            cid = self.get_drug_cid(drug_name)
            if cid:
                response['cid'] = cid
                molecular_weight = self.get_molecular_weight(cid)
                if molecular_weight:
                    response['molecular_weight'] = molecular_weight
            
            # Get hERG data from PubChem
            herg_ic50, herg_source = self.search_herg_data(drug_name)
            if herg_ic50 is not None:
                response['herg_ic50'] = herg_ic50
                response['herg_source'] = herg_source
            
            # Calculate concentrations if we have molecular weight
            if response['molecular_weight'] is not None:
                concentrations = self.calculate_concentrations(response['molecular_weight'], doses)
                
                # Format concentrations for output
                formatted_concentrations = []
                for conc in concentrations:
                    conc_data = {
                        'dose': conc.dose,
                        'theoretical_max': conc.theoretical_max,
                        'plasma_concentration': conc.plasma_concentration,
                        'ratio_theoretical': None,
                        'ratio_plasma': None
                    }
                    
                    # Calculate ratios if we have hERG IC50
                    if response['herg_ic50'] is not None:
                        try:
                            herg_ic50_float = float(str(response['herg_ic50']).strip())
                            if (conc.theoretical_max is not None and 
                                conc.plasma_concentration is not None):
                                if conc.theoretical_max > 0:
                                    conc_data['ratio_theoretical'] = round(
                                        herg_ic50_float / conc.theoretical_max, 3
                                    )
                                if conc.plasma_concentration > 0:
                                    conc_data['ratio_plasma'] = round(
                                        herg_ic50_float / conc.plasma_concentration, 3
                                    )
                        except (ValueError, TypeError, ZeroDivisionError) as e:
                            logger.warning(f"Error calculating ratios: {str(e)}")
                    
                    formatted_concentrations.append(conc_data)
                
                response['concentrations'] = formatted_concentrations
                
                # Determine theoretical binding
                if response['herg_ic50'] is not None:
                    response['theoretical_binding'] = any(
                        c['ratio_plasma'] is not None and c['ratio_plasma'] < 1 
                        for c in formatted_concentrations
                    )
            
            # Generate citations
            citations = []
            if response['cid']:
                citations.append(
                    "National Center for Biotechnology Information (2024). PubChem Compound Summary for CID {}, {}. "
                    "Retrieved January 10, 2024, from https://pubchem.ncbi.nlm.nih.gov/compound/{}".format(
                        response['cid'], drug_name, response['cid']
                    )
                )
            if response['herg_source']:
                citations.append(f"hERG IC50 data source: {response['herg_source']}")
            response['citations'] = citations
            
            return response
            
        except Exception as e:
            logger.error(f"Error analyzing drug: {str(e)}")
            return {
                'error': f"Error analyzing hERG activity: {str(e)}",
                'concentrations': []  # Ensure this is always present
            }
