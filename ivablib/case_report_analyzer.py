"""Module for analyzing case reports from medical papers."""

import re
import logging
import pandas as pd
from typing import List, Dict, Any, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CaseReportAnalyzer:
    """Analyzes medical case reports for relevant information."""
    
    def __init__(self):
        """Initialize the analyzer with regex patterns."""
        # Compile regex patterns for better performance
        self.patterns = {
            'age': re.compile(r'(\d+)[\s-]*(year|yr|y)[s\s-]*old|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
            'sex': re.compile(r'\b(male|female|man|woman|boy|girl)\b', re.IGNORECASE),
            'qtc': re.compile(r'QTc[\s:]*(\d+)', re.IGNORECASE),
            'qt_uncorrected': re.compile(r'\bQT[\s:]*(\d+)', re.IGNORECASE),
            'heart_rate': re.compile(r'(?:heart rate|HR|pulse)[\s:]*(\d+)(?:\s*beats?\s*per\s*min(?:ute)?|\s*bpm)?', re.IGNORECASE),
            'blood_pressure': re.compile(r'(?:blood pressure|BP)[\s:]*([\d/]+)', re.IGNORECASE),
            'had_tdp': re.compile(r'\b(?:torsade[s]* de pointes|TdP|torsades)\b', re.IGNORECASE),
            'patient_type': re.compile(r'\b(pediatric|child|infant|newborn|adult|elderly)\b', re.IGNORECASE),
            'outcome': re.compile(r'(?:treatment|therapy)[\s\w]+(?:successful|effective|improved|resolved|decreased|controlled|normalized)', re.IGNORECASE),
            'treatment_duration': re.compile(r'(?:after|within|for)\s*(\d+)\s*(?:hour|hr|h|day|week|month|year)s?', re.IGNORECASE),
            'drug_combination': re.compile(r'(?:with|plus|and|combination)\s+(amiodarone|flecainide|beta.?blocker|metoprolol|propranolol|sotalol)', re.IGNORECASE)
        }
        
    def extract_value(self, text: str, pattern: re.Pattern) -> Optional[str]:
        """Extract a value from text using a compiled regex pattern."""
        if not text:
            return None
            
        match = pattern.search(text)
        if match:
            # Get all groups that matched
            groups = [g for g in match.groups() if g is not None]
            if groups:
                # Return the first non-None group
                return groups[0].strip()
        return None
        
    def extract_all_matches(self, text: str, pattern: re.Pattern) -> List[str]:
        """Extract all matches from text using a compiled regex pattern."""
        if not text:
            return []
            
        matches = pattern.finditer(text)
        results = []
        for m in matches:
            groups = [g for g in m.groups() if g is not None]
            if groups:
                results.append(groups[0].strip())
        return results
        
    def analyze_paper(self, paper: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Analyze a single paper for case report information."""
        try:
            # Combine title and abstract for analysis
            text = f"{paper.get('title', '')} {paper.get('abstract', '')}"
            
            # Skip if no meaningful text to analyze
            if not text.strip():
                return None
            
            # Extract basic information
            age = self.extract_value(text, self.patterns['age'])
            sex = self.extract_value(text, self.patterns['sex'])
            qtc = self.extract_value(text, self.patterns['qtc'])
            qt = self.extract_value(text, self.patterns['qt_uncorrected'])
            hr = self.extract_value(text, self.patterns['heart_rate'])
            bp = self.extract_value(text, self.patterns['blood_pressure'])
            
            # Extract additional information
            patient_type = self.extract_value(text, self.patterns['patient_type'])
            outcome = bool(self.patterns['outcome'].search(text))
            treatment_duration = self.extract_value(text, self.patterns['treatment_duration'])
            drug_combinations = self.extract_all_matches(text, self.patterns['drug_combination'])
            
            # Check for TdP
            had_tdp = bool(self.patterns['had_tdp'].search(text))
            
            # Clean up numeric values
            try:
                age = float(age) if age and age.isdigit() else None
            except (ValueError, TypeError):
                age = None
                
            try:
                qtc = float(qtc) if qtc and qtc.isdigit() else None
            except (ValueError, TypeError):
                qtc = None
                
            try:
                qt = float(qt) if qt and qt.isdigit() else None
            except (ValueError, TypeError):
                qt = None
                
            try:
                hr = float(hr) if hr and hr.isdigit() else None
            except (ValueError, TypeError):
                hr = None
            
            # Only return if we found some relevant information
            if any([age, sex, qtc, qt, hr, bp, had_tdp, patient_type, outcome, treatment_duration, drug_combinations]):
                return {
                    'title': paper.get('title', ''),
                    'authors': paper.get('authors', ''),
                    'year': paper.get('year'),
                    'journal': paper.get('journal', ''),
                    'doi': paper.get('doi', ''),
                    'pmid': paper.get('pmid', ''),
                    'age': age,
                    'sex': sex.title() if sex else None,
                    'qtc': qtc,
                    'qt_uncorrected': qt,
                    'heart_rate': hr,
                    'blood_pressure': bp,
                    'had_tdp': 'Yes' if had_tdp else 'No',
                    'patient_type': patient_type.title() if patient_type else None,
                    'treatment_successful': 'Yes' if outcome else 'No',
                    'treatment_duration': treatment_duration,
                    'drug_combinations': ', '.join(drug_combinations) if drug_combinations else None
                }
            return None
            
        except Exception as e:
            logger.error(f"Error analyzing paper: {str(e)}")
            return None
            
    def analyze_papers(self, papers: List[Dict[str, Any]], drug_name: str) -> pd.DataFrame:
        """Analyze a list of papers and return results as a DataFrame."""
        try:
            results = []
            for paper in papers:
                result = self.analyze_paper(paper)
                if result:
                    results.append(result)
                    
            # Convert to DataFrame
            df = pd.DataFrame(results) if results else pd.DataFrame(columns=['title', 'authors', 'year', 'journal', 'doi', 'pmid', 'age', 'sex', 'qtc', 'qt_uncorrected', 'heart_rate', 'blood_pressure', 'had_tdp', 'patient_type', 'treatment_successful', 'treatment_duration', 'drug_combinations', 'drug'])
            
            # Add drug name
            if not df.empty:
                df['drug'] = drug_name
            
            return df
            
        except Exception as e:
            logger.error(f"Error analyzing papers: {str(e)}")
            return pd.DataFrame(columns=['title', 'authors', 'year', 'journal', 'doi', 'pmid', 'age', 'sex', 'qtc', 'qt_uncorrected', 'heart_rate', 'blood_pressure', 'had_tdp', 'patient_type', 'treatment_successful', 'treatment_duration', 'drug_combinations', 'drug'])

def analyze_papers(papers: List[Dict[str, Any]], drug_name: str) -> pd.DataFrame:
    """Analyze a list of papers and create a case report table."""
    analyzer = CaseReportAnalyzer()
    return analyzer.analyze_papers(papers, drug_name)

def convert_pdf_to_text(pdf_path: str) -> Optional[str]:
    """Convert a PDF file to text and save it"""
    try:
        # Create text file path
        text_path = pdf_path.rsplit('.', 1)[0] + '.txt'
        
        # Skip if text file already exists
        if os.path.exists(text_path):
            logger.info(f"Text file already exists: {text_path}")
            return text_path
            
        logger.info(f"Converting PDF to text: {pdf_path}")
        
        # Check if PDF exists
        if not os.path.exists(pdf_path):
            logger.error(f"PDF file not found: {pdf_path}")
            return None
            
        # Read PDF
        with open(pdf_path, 'rb') as file:
            try:
                reader = PyPDF2.PdfReader(file)
                if not reader.pages:
                    logger.error(f"PDF has no pages: {pdf_path}")
                    return None
                    
                text = ''
                page_count = 0
                
                # Extract text from each page
                for i, page in enumerate(reader.pages):
                    try:
                        page_text = page.extract_text()
                        if page_text:
                            # Clean and format the text
                            page_text = clean_text(page_text)
                            if page_text.strip():
                                text += f"\n[Page {i+1}]\n{page_text}\n"
                                page_count += 1
                                
                    except Exception as e:
                        logger.warning(f"Error extracting text from page {i+1}: {e}")
                        continue
                
                # Check if we got any text
                if not text.strip():
                    logger.error(f"No text could be extracted from any pages in {pdf_path}")
                    return None
                    
                logger.info(f"Successfully extracted text from {page_count} pages")
                
                # Save text to file
                with open(text_path, 'w', encoding='utf-8') as f:
                    f.write(text.strip())
                logger.info(f"Saved text to: {text_path}")
                return text_path
                    
            except Exception as e:
                logger.error(f"Error reading PDF {pdf_path}: {e}")
                return None
                
    except Exception as e:
        logger.error(f"Error converting PDF to text: {e}")
        return None
