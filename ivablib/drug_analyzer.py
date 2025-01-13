"""DrugAnalyzer class for analyzing drug-related literature."""

import logging
import pandas as pd
from typing import List, Dict, Any, Optional
import concurrent.futures
from .pubmed4125 import search_pubmed_case_reports, analyze_literature
from .case_report_analyzer import CaseReportAnalyzer

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DrugAnalyzer:
    """Analyzes drug-related literature and case reports."""
    
    def __init__(self):
        """Initialize the drug analyzer."""
        self.case_report_analyzer = CaseReportAnalyzer()
    
    def search_pubmed_case_reports(self, drug_name: str) -> pd.DataFrame:
        """Search PubMed for case reports about a drug."""
        return search_pubmed_case_reports(drug_name)
    
    def analyze_drug(self, drug_name: str) -> Dict[str, Any]:
        """Analyze literature for a given drug."""
        try:
            # Search for case reports
            case_reports = self.search_pubmed_case_reports(drug_name)
            if case_reports.empty:
                logger.warning(f"No case reports found for {drug_name}")
                return {}
                
            # Analyze the literature
            analysis = analyze_literature(drug_name)
            
            # Combine results
            result = {
                'drug_name': drug_name,
                'case_reports': case_reports.to_dict('records'),
                'analysis': analysis
            }
            
            return result
            
        except Exception as e:
            logger.error(f"Error analyzing drug: {str(e)}")
            return {}
