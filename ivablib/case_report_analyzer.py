"""
Case report analysis module for extracting structured information from medical case reports
using a language model.
"""

import os
import logging
import re
from typing import List, Dict, Optional
from dataclasses import dataclass, asdict
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM, pipeline
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set cache directory for models
os.environ['TRANSFORMERS_CACHE'] = '/opt/render/project/src/.cache/huggingface'
os.environ['HF_HOME'] = '/opt/render/project/src/.cache/huggingface'

@dataclass
class CaseReport:
    """Data class for storing case report information."""
    title: Optional[str] = None
    doi: Optional[str] = None 
    authors: Optional[str] = None
    age: Optional[str] = None
    sex: Optional[str] = None
    oral_dose: Optional[str] = None
    qt_uncorrected: Optional[str] = None
    qtc: Optional[str] = None
    heart_rate: Optional[str] = None
    had_tdp: Optional[str] = None
    blood_pressure: Optional[str] = None
    medical_history: Optional[str] = None
    medication_history: Optional[str] = None
    treatment_course: Optional[str] = None

class CaseReportAnalyzer:
    """Class for analyzing medical case reports using language models."""
    
    def __init__(self, model_name: str = "distilbert-base-cased-distilled-squad"):
        """Initialize the case report analyzer.
        
        Args:
            model_name: Name of the model to use
        """
        try:
            logger.info(f"Loading model {model_name}")
            self.qa_pipeline = pipeline(
                "question-answering",
                model=model_name,
                tokenizer=model_name,
                device=-1  # Force CPU
            )
            logger.info("Model loaded successfully")
        except Exception as e:
            logger.error(f"Error loading model: {str(e)}")
            raise

    def extract_answer(self, question: str, text: str) -> str:
        """Extract answer from text using question-answering pipeline.
        
        Args:
            question: Question to ask
            text: Text to extract answer from
            
        Returns:
            Answer string
        """
        try:
            result = self.qa_pipeline(
                question=question,
                context=text,
                max_length=128,
                truncation=True
            )
            return result['answer']
        except Exception as e:
            logger.error(f"Error extracting answer: {str(e)}")
            return None

    def clean_field(self, field: str, value: str) -> str:
        """Clean and validate field values.
        
        Args:
            field: Field name
            value: Value to clean
            
        Returns:
            Cleaned value
        """
        if not value:
            return None
            
        value = value.strip()
        
        # Field-specific cleaning
        if field == 'age':
            # Extract numeric age and standardize to years
            match = re.search(r'(\d+)\s*(year|yr|y|month|mo|m|week|wk|w|day|d)s?\s*old?', value.lower())
            if match:
                num, unit = match.groups()
                num = int(num)
                if unit.startswith(('month', 'mo', 'm')):
                    num = round(num/12, 1)
                elif unit.startswith(('week', 'wk', 'w')):
                    num = round(num/52, 1)
                elif unit.startswith(('day', 'd')):
                    num = round(num/365, 1)
                value = f"{num} years"
            
        elif field == 'sex':
            # Standardize sex/gender
            value = value.lower()
            if 'female' in value or 'woman' in value or 'f' in value:
                value = 'Female'
            elif 'male' in value or 'man' in value or 'm' in value:
                value = 'Male'
                
        elif field == 'oral_dose':
            # Extract numeric dose and unit
            match = re.search(r'(\d+\.?\d*)\s*(mg|mcg|g)', value.lower())
            if match:
                num, unit = match.groups()
                value = f"{num} {unit}"
                
        elif field == 'qt_uncorrected' or field == 'qtc':
            # Extract QT interval in milliseconds
            match = re.search(r'(\d+\.?\d*)\s*(ms|msec|millisecond)', value.lower())
            if match:
                num = match.group(1)
                value = f"{num} ms"
                
        elif field == 'heart_rate':
            # Extract heart rate in beats per minute
            match = re.search(r'(\d+\.?\d*)\s*(bpm|beat|b\.p\.m|per min)', value.lower())
            if match:
                num = match.group(1)
                value = f"{num} bpm"
                
        elif field == 'had_tdp':
            # Standardize TdP presence
            value = value.lower()
            if any(x in value for x in ['yes', 'had', 'developed', 'experienced', 'documented']):
                value = 'Yes'
            elif any(x in value for x in ['no', 'not', 'never', 'without']):
                value = 'No'
            else:
                value = 'Unknown'
                
        elif field == 'blood_pressure':
            # Extract systolic/diastolic in mmHg
            match = re.search(r'(\d+)/(\d+)\s*(mmhg|mm hg)?', value.lower())
            if match:
                sys, dia = match.groups()[:2]
                value = f"{sys}/{dia} mmHg"
                
        return value

    def analyze_case_report(self, text: str) -> CaseReport:
        """Analyze a case report and extract structured information.
        
        Args:
            text: Case report text
            
        Returns:
            CaseReport object with extracted information
        """
        # Create empty case report
        case = CaseReport()
        
        # Extract information for each field
        for field, prompt in self.prompts.items():
            try:
                # Get answer using QA pipeline
                answer = self.extract_answer(prompt, text)
                
                # Clean and validate the answer
                cleaned_value = self.clean_field(field, answer)
                
                # Set field value
                setattr(case, field, cleaned_value)
                
            except Exception as e:
                logger.error(f"Error extracting {field}: {str(e)}")
                setattr(case, field, None)
                
        return case

    def analyze_case_reports(self, texts: List[str]) -> pd.DataFrame:
        """Analyze multiple case reports and return results as DataFrame.
        
        Args:
            texts: List of case report texts
            
        Returns:
            DataFrame with extracted information
        """
        # Analyze each case report
        cases = []
        for text in texts:
            case = self.analyze_case_report(text)
            cases.append(asdict(case))
            
        # Convert to DataFrame
        df = pd.DataFrame(cases)
        
        return df

    def create_case_report_table(self, papers: List[object]) -> pd.DataFrame:
        """Create a table of case reports from papers.
        
        Args:
            papers: List of papers to analyze
            
        Returns:
            DataFrame with case report information
        """
        # Extract case info from papers
        all_cases = []
        for paper in papers:
            if paper.full_text:  # Access attribute directly
                case = self.analyze_case_report(paper.full_text)
                all_cases.append(asdict(case))
                
        # Convert to DataFrame
        if not all_cases:
            # Create empty DataFrame with all expected columns
            return pd.DataFrame(columns=[
                'title', 'doi', 'authors', 'age', 'sex', 'oral_dose', 'qt_uncorrected', 'qtc',
                'heart_rate', 'had_tdp', 'blood_pressure',
            ])
            
        # Create DataFrame from cases
        df = pd.DataFrame(all_cases)
        
        # Ensure all expected columns exist
        expected_cols = [
            'title', 'doi', 'authors', 'age', 'sex', 'oral_dose', 'qt_uncorrected', 'qtc',
            'heart_rate', 'had_tdp', 'blood_pressure',
        ]
        for col in expected_cols:
            if col not in df.columns:
                df[col] = None
                
        # Convert numeric columns
        numeric_cols = ['age', 'oral_dose', 'qt_uncorrected', 'qtc', 'heart_rate']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                
        return df

    def get_dose_patterns(self, drug_name: str) -> List[str]:
        """Get dose extraction patterns for a specific drug"""
        # Get drug names from CSV
        from .pubmed4125 import get_drug_names
        drug_names = get_drug_names(drug_name)
        if not drug_names:
            drug_names = [drug_name]
        drug_pattern = '|'.join(drug_names)  # e.g. "galantamine|razadyne|reminyl"
        
        return [
            f'(?:{drug_pattern}).*?(\\d+(?:\\.\\d+)?)\s*mg',
            f'(\\d+(?:\\.\\d+)?)\s*mg.*?(?:{drug_pattern})',
            r'(?:dose|dosage|prescribed).*?(\d+(?:\.\d+)?)\s*mg',
            r'(\d+(?:\.\d+)?)\s*mg\s*(?:bid|b\.i\.d|twice daily)',
            r'treated with.*?(\d+(?:\.\d+)?)\s*mg'
        ]

    def analyze_papers(self, papers: List[Dict], drug_name: Optional[str] = None) -> pd.DataFrame:
        """Analyze a list of papers and return a pandas DataFrame with the results.
        
        Args:
            papers: List of paper dictionaries
            drug_name: Name of the drug being analyzed
            
        Returns:
            DataFrame with case report information
        """
        if not papers:
            logger.warning("No papers to analyze")
            return pd.DataFrame(columns=[
                'title', 'doi', 'pmid', 'had_tdp', 'age', 'sex', 'oral_dose',
                'heart_rate', 'blood_pressure', 'qtc', 'max_qtc', 'baseline_qtc'
            ])
        
        case_info = []
        
        # Get dose patterns for this drug
        dose_patterns = self.get_dose_patterns(drug_name) if drug_name else []
        
        for paper in papers:
            try:
                # Initialize case with basic info
                case = {
                    'title': paper.get('TI', ''),
                    'doi': paper.get('DOI', ''),  # Use the DOI we extracted earlier
                    'pmid': paper.get('PMID', ''),
                    'drug_name': drug_name,
                    'text_source': 'full_text' if paper.get('full_text') else 'abstract',
                    'had_tdp': False,
                    'age': None,
                    'sex': '',
                    'oral_dose': None,
                    'heart_rate': None,
                    'blood_pressure': '',
                    'qtc': None,
                    'max_qtc': None,
                    'baseline_qtc': None
                }
                
                # Get text content - prefer full text over abstract
                text = paper.get('full_text', '') or paper.get('AB', '')
                if not text:
                    logger.warning(f"No text content found for paper {case['title']}")
                    continue
                
                text_lower = text.lower()
                
                # Extract oral dose
                for pattern in dose_patterns:
                    dose_match = re.search(pattern, text_lower)
                    if dose_match:
                        try:
                            case['oral_dose'] = float(dose_match.group(1))
                            break
                        except (ValueError, IndexError):
                            continue
                
                # Extract age with simple, reliable patterns
                age_patterns = [
                    # Direct age mentions with word boundaries
                    r'(?<!\d)(?:aged?|age[d\s]of)\s*(\d{2,3})(?!\d)',
                    r'(?<!\d)(\d{2,3})[-\s]*(?:years?|yrs?|y\.?o\.?)[-\s]*old(?!\d)',
                    r'age\s*(?:was|is|:)\s*(\d{2,3})(?!\d)',
                    # Age with patient descriptors
                    r'(?:patient|man|woman|male|female|person|case).*?(?:aged?|age)\s*(\d{2,3})(?!\d)',
                    r'(?:patient|man|woman|male|female|person|case).*?(\d{2,3})\s*(?:years?|yrs?|y\.?o\.?)(?!\d)',
                    # Age in parentheses
                    r'\((?:aged?|age)?\s*(\d{2,3})\s*(?:years?|yrs?|y\.?o\.?)?\)',
                    # Age at start of sentence
                    r'(?:^|\.)\s*(?:a|an?)\s*(\d{2,3})[-\s]*(?:years?|yrs?|y\.?o\.?)[-\s]*old'
                ]
                
                # Search in both title and text
                text_to_search = f"{text_lower} {case['title'].lower()}"
                logger.debug(f"Searching for age in text: {text_to_search[:200]}...")  # Show first 200 chars
                
                found_ages = []
                for pattern in age_patterns:
                    for age_match in re.finditer(pattern, text_to_search):
                        try:
                            age = int(age_match.group(1))
                            if 10 <= age <= 100:  # Basic validation
                                found_ages.append(age)
                                context = text_to_search[max(0, age_match.start()-30):min(len(text_to_search), age_match.end()+30)]
                                logger.debug(f"Found potential age: {age} in context: ...{context}...")
                        except ValueError:
                            continue
                
                # If we found ages, prefer the one in elderly range or most common
                if found_ages:
                    # First try to find ages in the elderly range (65-95)
                    elderly_ages = [age for age in found_ages if 65 <= age <= 95]
                    if elderly_ages:
                        # If multiple elderly ages, use the most common one
                        age_counts = {}
                        for age in elderly_ages:
                            age_counts[age] = elderly_ages.count(age)
                        case['age'] = max(age_counts.items(), key=lambda x: x[1])[0]
                    else:
                        # If no elderly ages, use the most common age
                        age_counts = {}
                        for age in found_ages:
                            age_counts[age] = found_ages.count(age)
                        case['age'] = max(age_counts.items(), key=lambda x: x[1])[0]
                    
                    logger.info(f"Selected age {case['age']} from candidates: {found_ages}")
                else:
                    logger.debug(f"No age found in text: {text_to_search[:200]}...")
                
                # Extract sex with improved patterns
                sex_patterns = {
                    'woman': [
                        r'\b(?:female|woman|lady)\b',
                        r'\bshe\b.*?patient',
                        r'woman.*?patient',
                        r'female.*?patient',
                        r'lady.*?patient'
                    ],
                    'man': [
                        r'\b(?:male|man|gentleman)\b',
                        r'\bhe\b.*?patient',
                        r'man.*?patient',
                        r'male.*?patient',
                        r'gentleman.*?patient'
                    ]
                }
                
                # Look for sex indicators in title and text
                text_to_search = text_lower + ' ' + case['title'].lower()
                for sex, patterns in sex_patterns.items():
                    if any(re.search(pattern, text_to_search) for pattern in patterns):
                        case['sex'] = sex
                        logger.info(f"Found sex: {sex}")
                        break
                
                # Extract heart rate with improved patterns
                hr_patterns = [
                    r'(?:heart rate|hr|pulse|rate)\s*(?:was|of|[:=]|\s)\s*(\d{2,3})\s*(?:bpm|beats|b\.?p\.?m\.?)?',
                    r'(?:heart rate|hr|pulse|rate)\s*(?:of|[:=]|\s)\s*(\d{2,3})\s*(?:bpm|beats|b\.?p\.?m\.?)?',
                    r'(?:bradycardia|tachycardia).*?(\d{2,3})\s*(?:bpm|beats|b\.?p\.?m\.?)?',
                    r'(?:pulse|rate)\s*(\d{2,3})\s*(?:bpm|beats|b\.?p\.?m\.?)?',
                    r'(?:heart\s+rate|hr)\s*(?:increased|decreased)\s*to\s*(\d{2,3})',
                    r'(?:heart\s+rate|hr)\s*(?:of|at)\s*(\d{2,3})'
                ]
                
                for pattern in hr_patterns:
                    hr_match = re.search(pattern, text_lower)
                    if hr_match:
                        try:
                            hr = int(hr_match.group(1))
                            if 20 <= hr <= 200:  # Basic validation
                                case['heart_rate'] = hr
                                logger.info(f"Found heart rate: {hr}")
                                break
                        except ValueError:
                            continue
            
                # Extract blood pressure with improved patterns
                bp_patterns = [
                    r'(?:blood pressure|bp)\s*(?:was|of|[:=]|\s)\s*(\d{2,3})[/\\](\d{2,3})',
                    r'(?:blood pressure|bp)\s*(?:of|[:=]|\s)\s*(\d{2,3})[/\\](\d{2,3})',
                    r'(\d{2,3})[/\\](\d{2,3})\s*(?:mmhg|mm\s*hg)',
                    r'(?:systolic|diastolic).*?(\d{2,3})[/\\](\d{2,3})',
                    r'bp\s*(?:was|:)?\s*(\d{2,3})[/\\](\d{2,3})',
                    r'(?:pressure|bp)\s*(?:of|at)\s*(\d{2,3})[/\\](\d{2,3})'
                ]
                
                for pattern in bp_patterns:
                    bp_match = re.search(pattern, text_lower)
                    if bp_match:
                        try:
                            systolic = int(bp_match.group(1))
                            diastolic = int(bp_match.group(2))
                            if 60 <= systolic <= 220 and 40 <= diastolic <= 120:  # Basic validation
                                case['blood_pressure'] = f"{systolic}/{diastolic}"
                                logger.info(f"Found blood pressure: {systolic}/{diastolic}")
                                break
                        except ValueError:
                            continue
                
                # Extract QTc values
                qtc_patterns = [
                    r'qtc.*?(\d{3})(?:\s*ms)?',
                    r'qt[c]?\s*(?:interval|prolongation).*?(\d{3})(?:\s*ms)?',
                    r'(?:bazett|fridericia).*?(\d{3})(?:\s*ms)?',
                    r'(?:prolonged|increased)\s*qtc.*?(\d{3})(?:\s*ms)?',
                    r'(?:qtc|qt)\s*(?:of|was|=|:|\s)\s*(\d{3})(?:\s*ms)?',
                    r'(?:qtc|qt)\s*(?:interval)?\s*(?:of|was|=|:|\s)\s*(\d{3})(?:\s*ms)?'
                ]
                
                qtc_values = []
                for pattern in qtc_patterns:
                    for match in re.finditer(pattern, text_lower):
                        try:
                            qtc = int(match.group(1))
                            if 300 <= qtc <= 700:  # Basic validation
                                qtc_values.append(qtc)
                                logger.info(f"Found QTc value: {qtc} ms")
                        except ValueError:
                            continue
                
                if qtc_values:
                    case['qtc'] = qtc_values[0]  # First QTc mentioned
                    case['max_qtc'] = max(qtc_values)
                    case['baseline_qtc'] = min(qtc_values)
                    logger.info(f"QTc values found - First: {case['qtc']}, Max: {case['max_qtc']}, Baseline: {case['baseline_qtc']}")
                
                # Check for TdP
                tdp_patterns = [
                    r'\b(?:torsade|torsades|tdp)\b',
                    r'polymorphic\s+ventricular\s+tachycardia',
                    r'(?:developed|experienced|had).*?torsade',
                    r'torsade.*?(?:occurred|developed|observed)',
                    r'(?:episodes?|events?)\s+of\s+(?:torsade|torsades|tdp)'
                ]
                
                case['had_tdp'] = any(re.search(pattern, text_lower) for pattern in tdp_patterns)
                
                case_info.append(case)
            
            except Exception as e:
                logger.error(f"Error analyzing paper: {str(e)}")
                continue
    
        # Create DataFrame
        if case_info:
            df = pd.DataFrame(case_info)
            
            # Save to Excel if drug name provided
            if drug_name:
                # Create output directory if it doesn't exist
                output_dir = os.path.join(os.getcwd(), 'papers', drug_name.lower())
                os.makedirs(output_dir, exist_ok=True)
                
                # Get current time in format HHMM
                from datetime import datetime
                timestamp = datetime.now().strftime("%H%M")
                
                # Save Excel file with timestamp
                output_file = os.path.join(output_dir, f"{drug_name.lower()}_{timestamp}.xlsx")
                df.to_excel(output_file, index=False)
                logger.info(f"Saved analysis to {output_file}")
            
            return df
        else:
            logger.warning("No papers could be analyzed")
            return pd.DataFrame()

def analyze_papers(papers: List[object], drug_name: Optional[str] = None) -> pd.DataFrame:
    """Analyze a list of papers and create a case report table.
    
    Args:
        papers: List of papers to analyze
        drug_name: Name of the drug being analyzed
        
    Returns:
        DataFrame with case report information
    """
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
