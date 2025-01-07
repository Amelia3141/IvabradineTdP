"""Drug concentration and risk calculations."""
from typing import Dict, Optional, Tuple
from .chembl import fetch_pubchem_data, get_herg_ic50

def calculate_concentrations(dose_mg: float, molecular_weight: float = 427.5,
                           volume_of_distribution: float = 252,
                           bioavailability: float = 0.4) -> Dict[str, float]:
    """Calculate theoretical and bioavailable concentrations."""
    # Convert dose to moles
    dose_mol = dose_mg / molecular_weight / 1000  # Convert to moles
    
    # Calculate theoretical concentration
    theoretical_conc = dose_mol * 1e6 / volume_of_distribution  # Convert to μM
    
    # Calculate bioavailable concentration
    bioavailable_conc = theoretical_conc * bioavailability
    
    return {
        'theoretical': round(theoretical_conc, 3),
        'bioavailable': round(bioavailable_conc, 3)
    }

def calculate_drug_concentrations(drug_name: str, dose_mg: float) -> Dict[str, Optional[float]]:
    """Calculate drug concentrations and hERG ratios."""
    results = {
        'theoretical_conc': None,
        'bioavailable_conc': None,
        'herg_ic50': None,
        'safety_margin': None
    }
    
    # Get molecular properties
    pubchem_data = fetch_pubchem_data(drug_name)
    if not pubchem_data or 'molecular_weight' not in pubchem_data:
        return results
        
    # Get hERG IC50
    herg_ic50 = get_herg_ic50(drug_name)
    if not herg_ic50:
        return results
        
    # Calculate concentrations
    concs = calculate_concentrations(
        dose_mg=dose_mg,
        molecular_weight=pubchem_data['molecular_weight']
    )
    
    results['theoretical_conc'] = concs['theoretical']
    results['bioavailable_conc'] = concs['bioavailable']
    results['herg_ic50'] = herg_ic50
    
    # Calculate safety margin
    if concs['bioavailable'] > 0:
        results['safety_margin'] = herg_ic50 / concs['bioavailable']
    
    return results
