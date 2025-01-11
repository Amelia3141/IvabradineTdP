"""Module for analyzing case reports from medical papers."""

import re
import logging
import pandas as pd
from typing import List, Dict, Any, Optional
import os
import PyPDF2

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CaseReportAnalyzer:
    """Analyzes medical case reports for relevant information."""
    
    def __init__(self):
        """Initialize with regex patterns for extracting information."""
        self.patterns = {
            'age': re.compile(r'(?:age[d]?|year[s]?[- ]old)[: ]*(\d+)|\b(\d+)[ -](?:year[s]?[- ]old|yo)\b', re.IGNORECASE),
            'sex': re.compile(r'\b(?:male|female|[MF])\b', re.IGNORECASE),
            'qtc': re.compile(r'QTc[: ]*(\d+)(?:\s*m?sec)?|corrected[- ]QT[: ]*(\d+)', re.IGNORECASE),
            'qt_uncorrected': re.compile(r'QT[: ]*(\d+)(?:\s*m?sec)?(?![Cc])', re.IGNORECASE),
            'heart_rate': re.compile(r'(?:heart rate|HR|pulse)[: ]*(\d+)(?:\s*(?:bpm|beats per minute|/min))?', re.IGNORECASE),
            'blood_pressure': re.compile(r'(?:blood pressure|BP)[: ]*(\d+/\d+)', re.IGNORECASE),
            'had_tdp': re.compile(r'\b(?:torsade[s]? de pointes|tdp|pmvt)\b(?! risk| possible| unlikely| negative| no | not | without)', re.IGNORECASE),
            'patient_type': re.compile(r'(?:inpatient|outpatient|emergency|hospitalized)', re.IGNORECASE),
            'treatment_duration': re.compile(r'(?:treated for|duration of|over|for)\s+(\d+)\s+(?:day|week|month|year)s?', re.IGNORECASE),
            'outcome': re.compile(r'(?:recovered|improved|resolved|discharged|survived|died|expired|fatal)', re.IGNORECASE),
            'drug_combination': re.compile(r'(?:with|plus|and|combination)\s+(amiodarone|flecainide|beta.?blocker|metoprolol|propranolol|sotalol)', re.IGNORECASE)
        }
        
    def _extract_heart_rate(self, text):
        """Extract heart rate from text."""
        hr_patterns = [
            r'(?:heart rate|HR|pulse)[: ]*(\d+)(?:\s*(?:bpm|beats per minute|/min))?',
            r'HR[: ]*(\d+)',
            r'pulse[: ]*(\d+)',
            r'(?:heart rate|HR)[- ](?:was|of)[: ]*(\d+)',
            r'(?:heart rate|HR)[- ](?:increased|decreased)[- ]to[: ]*(\d+)'
        ]
        
        for pattern in hr_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                try:
                    hr = int(match.group(1))
                    if 20 <= hr <= 300:  # Reasonable heart rate range
                        return hr
                except ValueError:
                    continue
        return None

    def _extract_qtc(self, text):
        """Extract QTc interval from text."""
        qtc_patterns = [
            r'QTc[: ]*(\d+)(?:\s*m?sec)?',
            r'corrected[- ]QT[: ]*(\d+)',
            r'QTc interval[: ]*(\d+)',
            r'QTc prolongation to[: ]*(\d+)',
            r'QTc was[: ]*(\d+)',
            r'prolonged QTc of[: ]*(\d+)'
        ]
        
        for pattern in qtc_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                try:
                    qtc = int(match.group(1))
                    if 300 <= qtc <= 700:  # Reasonable QTc range
                        return qtc
                except ValueError:
                    continue
        return None

    def _extract_age(self, text):
        """Extract age from text."""
        age_patterns = [
            r'(?:age[d]?|year[s]?[- ]old)[: ]*(\d+)',
            r'\b(\d+)[ -](?:year[s]?[- ]old|yo)\b',
            r'(\d+)[ -]year[- ]old',
            r'age[d]?:?\s*(\d+)',
            r'(\d+)(?:\s*(?:y|yr|year)s?\s*(?:old)?)'
        ]
        
        for pattern in age_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                try:
                    age = int(match.group(1))
                    if 0 <= age <= 120:  # Reasonable age range
                        return age
                except ValueError:
                    continue
        return None

    def _extract_sex(self, text):
        """Extract sex from text."""
        # Look for explicit mentions first
        if re.search(r'\b(?:male|man)\b', text, re.IGNORECASE):
            return 'Male'
        if re.search(r'\b(?:female|woman)\b', text, re.IGNORECASE):
            return 'Female'
            
        # Then look for single letter indicators
        if re.search(r'\bM\b', text):
            return 'Male'
        if re.search(r'\bF\b', text):
            return 'Female'
            
        return None

    def _extract_tdp_status(self, text):
        """Extract TdP status from text."""
        text = text.lower()
        
        # Check for explicit negatives first
        negative_patterns = [
            r'no[t]?\s+(?:evidence\s+of\s+)?(?:torsade|tdp)',
            r'without\s+(?:evidence\s+of\s+)?(?:torsade|tdp)',
            r'denied\s+(?:evidence\s+of\s+)?(?:torsade|tdp)',
            r'negative\s+for\s+(?:torsade|tdp)'
        ]
        
        for pattern in negative_patterns:
            if re.search(pattern, text):
                return 'Negative'
        
        # Then check for positives
        positive_patterns = [
            r'\b(?:torsade[s]?\s+de\s+pointes|tdp)\b',
            r'polymorphic\s+[vV]entricular\s+[tT]achycardia',
            r'\bpmvt\b',
            r'developed\s+(?:torsade|tdp)',
            r'experienced\s+(?:torsade|tdp)'
        ]
        
        for pattern in positive_patterns:
            if re.search(pattern, text):
                return 'Positive'
        
        return 'Not reported'

    def analyze_paper(self, paper):
        """Analyze a single paper for case report information."""
        try:
            # Try to get full text if available
            full_text = None
            if paper.get('pmid'):
                try:
                    from .pubmed4125 import get_full_text
                    full_text = get_full_text(paper['pmid'])
                except Exception as e:
                    logger.warning(f"Could not get full text for PMID {paper['pmid']}: {e}")
            
            # Use full text if available, otherwise use title + abstract
            text = full_text if full_text else f"{paper.get('title', '')} {paper.get('abstract', '')}"
            
            # Extract information
            age = self._extract_age(text)
            sex = self._extract_sex(text)
            heart_rate = self._extract_heart_rate(text)
            qtc = self._extract_qtc(text)
            tdp_status = self._extract_tdp_status(text)
            
            # Log what we found
            logger.info(f"Paper {paper.get('pmid')}: Age={age}, Sex={sex}, HR={heart_rate}, QTc={qtc}, TdP={tdp_status}")
            
            return {
                'pmid': paper.get('pmid'),
                'age': age,
                'sex': sex,
                'heart_rate': heart_rate,
                'qtc': qtc,
                'tdp_status': tdp_status
            }
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
            df = pd.DataFrame(results) if results else pd.DataFrame(columns=['pmid', 'age', 'sex', 'heart_rate', 'qtc', 'tdp_status'])
            
            # Add drug name
            if not df.empty:
                df['drug'] = drug_name
            
            return df
            
        except Exception as e:
            logger.error(f"Error analyzing papers: {str(e)}")
            return pd.DataFrame(columns=['pmid', 'age', 'sex', 'heart_rate', 'qtc', 'tdp_status', 'drug'])

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
