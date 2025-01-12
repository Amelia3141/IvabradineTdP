"""Module for analyzing case reports from medical papers."""

import re
import logging
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
import os
import PyPDF2

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def clean_text(text: str) -> str:
    """Clean and normalize text."""
    if not text:
        return ""
    
    # Remove excessive whitespace
    text = ' '.join(text.split())
    
    # Fix common OCR/PDF extraction issues
    text = text.replace('|', 'I')  # Common OCR error
    text = text.replace('—', '-')  # Normalize dashes
    text = text.replace('–', '-')
    text = text.replace(''', "'")
    text = text.replace('"', '"')
    text = text.replace('"', '"')
    
    # Fix spacing around numbers and units
    text = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', text)  # Add space between number and unit
    text = re.sub(r'(\d)\s+([.,])\s*(\d)', r'\1\2\3', text)  # Fix decimal numbers
    
    return text

class CaseReportAnalyzer:
    """Analyzes medical case reports for relevant information."""
    
    def __init__(self):
        """Initialize the analyzer with regex patterns."""
        # Compile regex patterns for better performance
        self.patterns = {
            # Enhanced demographic patterns with more variations
            'age': re.compile(r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
            
            'sex': re.compile(r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b', re.IGNORECASE),
            
            # Enhanced dosing patterns with medical abbreviations
            'oral_dose': re.compile(r'(?:oral|po|by\s+mouth|orally|per\s+os)\s*(?:dose|dosage|doses|administration|administered|given|taken|prescribed|started|initiated|received|treated\s+with)?\s*(?:with|at|of|a|total|daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid)?\s*(?:dose\s+of\s+)?(\d+\.?\d*)\s*(?:mg|milligrams?|g|grams?|mcg|micrograms?|µg)', re.IGNORECASE),
            
            # Enhanced ECG patterns with more variations
            'qtc': re.compile(r'(?:QTc[FBR]?|corrected\s+QT|QT\s*corrected|QT[c]?\s*interval\s*(?:corrected)?|corrected\s*QT\s*interval)[\s:]*(?:interval|duration|measurement|prolongation|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
            
            'qt': re.compile(r'\b(?:QT|uncorrected\s+QT|QT\s*interval|interval\s*QT|baseline\s*QT)[\s:]*(?:interval|duration|measurement|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
            
            'heart_rate': re.compile(r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:|decreased\s+to|increased\s+to)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?', re.IGNORECASE),
            
            'blood_pressure': re.compile(r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:|measured)?\s*([\d/]+)(?:\s*mm\s*Hg)?', re.IGNORECASE),
            
            # Enhanced arrhythmia patterns
            'tdp': re.compile(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', re.IGNORECASE),
            
            # Enhanced history patterns with better context capture
            'medical_history': re.compile(r'(?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            'medication_history': re.compile(r'(?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            # Enhanced treatment patterns
            'treatment_course': re.compile(r'(?:treatment|therapy|management|intervention|administered|given)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
        }

    def get_smart_context(self, text: str, match_start: int, match_end: int, min_chars: int = 50) -> str:
        """
        Get context with a minimum of min_chars but expanded to complete sentences.
        
        Args:
            text: The full text to extract context from
            match_start: Start position of the match
            match_end: End position of the match
            min_chars: Minimum number of characters for context window
            
        Returns:
            String containing the smart context window
        """
        # Get initial context window
        start_pos = max(0, match_start - min_chars)
        end_pos = min(len(text), match_end + min_chars)
        
        # Expand to nearest sentence boundaries
        # Look backwards for sentence start
        while start_pos > 0:
            prev_char = text[start_pos - 1]
            if prev_char in '.!?' and text[start_pos].isupper():
                break
            start_pos -= 1
        
        # Look forwards for sentence end
        while end_pos < len(text):
            if text[end_pos - 1] in '.!?' and (end_pos == len(text) or text[end_pos].isupper()):
                break
            end_pos += 1
        
        return text[start_pos:end_pos].strip()

    def extract_value_with_context(self, text: str, pattern: re.Pattern) -> Optional[Tuple[str, str]]:
        """Extract a value and its context from text using a compiled regex pattern."""
        if not text:
            return None
            
        match = pattern.search(text)
        if match:
            # Return the first non-None group
            return next((g for g in match.groups() if g is not None), None), self.get_smart_context(text, match.start(), match.end())
        return None

    def analyze_qt_comprehensive(self, text: str) -> Dict[str, Any]:
        """Perform comprehensive QT analysis on text with smart context windows."""
        results = {
            'numeric_values': {},
            'tdp_mentions': [],
            'drug_qt_interactions': [],
            'brady_qt_relations': [],
            'figure_qt': [],
            'table_qt': [],
            'supp_qt': []
        }
        
        # Find numeric values with units
        for pattern_name in ['qtc', 'qt']:
            matches = self.patterns[pattern_name].finditer(text)
            for match in matches:
                if pattern_name not in results['numeric_values']:
                    results['numeric_values'][pattern_name] = []
                value = float(match.group(1))
                # Convert seconds to milliseconds if needed
                if any(unit in match.group(0).lower() for unit in ['sec', 's', 'seconds']):
                    value *= 1000
                results['numeric_values'][pattern_name].append({
                    'value': value,
                    'context': self.get_smart_context(text, match.start(), match.end())
                })
        
        # Find Torsades mentions
        tdp_matches = self.patterns['tdp'].finditer(text)
        for match in tdp_matches:
            results['tdp_mentions'].append({
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        # Find drug-QT relationships
        drug_matches = self.patterns['drug_qt'].finditer(text)
        for match in drug_matches:
            results['drug_qt_interactions'].append({
                'drug': match.group('drug'),
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        # Find bradycardia-QT relationships
        brady_matches = self.patterns['brady_qt'].finditer(text)
        for match in brady_matches:
            results['brady_qt_relations'].append({
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        # Find figure/table QT mentions
        figure_matches = self.patterns['figure_qt'].finditer(text)
        for match in figure_matches:
            results['figure_qt'].append({
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        # Find table cell QT values
        table_matches = self.patterns['table_qt'].finditer(text)
        for match in table_matches:
            results['table_qt'].append({
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        # Find supplementary material QT values
        supp_matches = self.patterns['supp_qt'].finditer(text)
        for match in supp_matches:
            results['supp_qt'].append({
                'text': match.group(0).strip(),
                'context': self.get_smart_context(text, match.start(), match.end())
            })
        
        return results

    def analyze_papers(self, papers: List[Dict[str, Any]], drug_name: str) -> pd.DataFrame:
        """Analyze a list of papers and create a case report table."""
        # Define columns for the output DataFrame
        columns = [
            'Case Report Title',
            'Age', 
            'Sex',
            'Oral Dose (mg)',
            'theoretical max concentration (μM)',
            '40% bioavailability',
            'Theoretical HERG IC50 / Concentration μM',
            '40% Plasma concentration',
            'Uncorrected QT (ms)',
            'QTc',
            'QTR',
            'QTF',
            'Heart Rate (bpm)',
            'Torsades de Pointes?',
            'Blood Pressure (mmHg)',
            'Medical History',
            'Medication History',
            'Course of Treatment'
        ]
        
        results = []
        for paper in papers:
            try:
                # Extract text from paper
                text = paper.get('FullText', '') or paper.get('Abstract', '')
                if not text or len(text.strip()) < 100:  # Skip if text is too short
                    logger.warning(f"Skipping paper with insufficient text: {paper.get('Title', 'Untitled')}")
                    continue
                
                # Clean the text
                text = clean_text(text)
                
                # Extract all values with context
                result = {
                    'title': paper.get('Title', 'Untitled Paper'),
                    'pmid': paper.get('PMID', 'N/A'),
                    'doi': paper.get('DOI', ''),
                    'age': self._extract_with_context(text, self.patterns['age']),
                    'sex': self._clean_sex(self._extract_with_context(text, self.patterns['sex'])),
                    'oral_dose_value': self._extract_numeric(self._extract_with_context(text, self.patterns['oral_dose'])),
                    'oral_dose_unit': self._extract_unit(self._extract_with_context(text, self.patterns['oral_dose'])),
                    'oral_dose_freq': self._extract_frequency(self._extract_with_context(text, self.patterns['oral_dose'])),
                    'qt_value': self._extract_numeric(self._extract_with_context(text, self.patterns['qt'])),
                    'qtc_value': self._extract_numeric(self._extract_with_context(text, self.patterns['qtc'])),
                    'heart_rate_value': self._extract_numeric(self._extract_with_context(text, self.patterns['heart_rate'])),
                    'blood_pressure_value': self._extract_with_context(text, self.patterns['blood_pressure']),
                    'tdp_present': bool(self.patterns['tdp'].search(text)),
                    'tdp_context': self.get_smart_context(text, self.patterns['tdp'].search(text).start(), self.patterns['tdp'].search(text).end()) if self.patterns['tdp'].search(text) else None,
                    'medical_history': self._extract_with_context(text, self.patterns['medical_history']),
                    'medication_history': self._extract_with_context(text, self.patterns['medication_history']),
                    'treatment_course': self._extract_with_context(text, self.patterns['treatment_course'])
                }
                
                # Calculate QTR and QTF if possible
                qt_val = self.extract_numeric(result['qt_value'])
                hr_val = self.extract_numeric(result['heart_rate_value'])
                if qt_val and hr_val:
                    try:
                        # Calculate QTR using Rautaharju formula
                        qtr = qt_val / (120 / hr_val) ** 0.5
                        result['qtr'] = f"{qtr:.1f}"
                        
                        # Calculate QTF using Fridericia formula
                        qtf = qt_val / (60 / hr_val) ** (1/3)
                        result['qtf'] = f"{qtf:.1f}"
                    except Exception as e:
                        logger.error(f"Error calculating QTR/QTF: {e}")
                        result['qtr'] = ''
                        result['qtf'] = ''
                else:
                    result['qtr'] = ''
                    result['qtf'] = ''
                
                # Format result to match desired columns
                formatted_result = {
                    'Case Report Title': result['title'],
                    'Age': result['age'] or '',
                    'Sex': result['sex'] or '',
                    'Oral Dose (mg)': result['oral_dose_value'] or '',
                    'theoretical max concentration (μM)': result['theoretical_max'] or '',
                    '40% bioavailability': result['bioavailability'] or '',
                    'Theoretical HERG IC50 / Concentration μM': result['herg_ic50'] or '',
                    '40% Plasma concentration': result['plasma_concentration'] or '',
                    'Uncorrected QT (ms)': result['qt_value'] or '',
                    'QTc': result['qtc_value'] or '',
                    'QTR': result['qtr'] or '',
                    'QTF': result['qtf'] or '',
                    'Heart Rate (bpm)': result['heart_rate_value'] or '',
                    'Torsades de Pointes?': result['tdp_present'],
                    'Blood Pressure (mmHg)': result['blood_pressure_value'] or '',
                    'Medical History': result['medical_history'] or '',
                    'Medication History': result['medication_history'] or '',
                    'Course of Treatment': result['treatment_course'] or ''
                }
                
                results.append(formatted_result)
                logger.info(f"Successfully analyzed paper {paper.get('PMID', 'Unknown PMID')}")
                
            except Exception as e:
                logger.error(f"Error analyzing paper {paper.get('PMID', 'Unknown PMID')}: {e}")
                continue
        
        # Create DataFrame with specific columns
        df = pd.DataFrame(results, columns=columns)
        
        # Clean up the DataFrame
        df = df.replace({None: '', 'None': '', 'nan': ''})
        df = df.fillna('')
        
        return df
        
    def analyze_paper(self, paper: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Analyze a single paper for case report information."""
        try:
            # Combine full text and abstract
            text = paper.get('FullText', '') + ' ' + paper.get('Abstract', '')
            if not text.strip():
                return {}

            # Initialize report
            report = {
                'Title': paper.get('Title', ''),
                'PMID': paper.get('PMID', ''),
                'DOI': paper.get('DOI', ''),
                'Abstract': paper.get('Abstract', '')
            }

            # Define patterns for information extraction
            patterns = {
                'Age': re.compile(r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
                'Sex': re.compile(r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b', re.IGNORECASE),
                'Oral Dose (mg)': re.compile(r'(?:oral|po|by\s+mouth|orally|per\s+os)\s*(?:dose|dosage|doses|administration|administered|given|taken|prescribed|started|initiated|received|treated\s+with)?\s*(?:with|at|of|a|total|daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid)?\s*(?:dose\s+of\s+)?(\d+\.?\d*)\s*(?:mg|milligrams?|g|grams?|mcg|micrograms?|µg)', re.IGNORECASE),
                'Theoretical Max Concentration (μM)': re.compile(r'(?:theoretical|max|maximum|peak)\s*(?:concentration|level|plasma|blood)\s*(?:of|was|is|=|:)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar)', re.IGNORECASE),
                '40% Bioavailability': re.compile(r'(?:bioavailability|F|BA)\s*(?:of|was|is|=|:)?\s*(\d+\.?\d*)\s*%', re.IGNORECASE),
                'Theoretical HERG IC50 / Concentration μM': re.compile(r'(?:herg|ikr)\s*(?:ic50|ic\s*50)\s*(?:of|was|is|=|:)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar)', re.IGNORECASE),
                '40% Plasma Concentration': re.compile(r'(?:plasma|blood|serum)\s*(?:concentration|level)\s*(?:of|was|is|=|:)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar)', re.IGNORECASE),
                'Uncorrected QT (ms)': re.compile(r'\b(?:QT|uncorrected\s+QT|QT\s*interval|interval\s*QT|baseline\s*QT)[\s:]*(?:interval|duration|measurement|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
                'QTc': re.compile(r'(?:QTc[FBR]?|corrected\s+QT|QT\s*corrected|QT[c]?\s*interval\s*(?:corrected)?|corrected\s*QT\s*interval)[\s:]*(?:interval|duration|measurement|prolongation|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
                'Heart Rate (bpm)': re.compile(r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:|decreased\s+to|increased\s+to)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?', re.IGNORECASE),
                'Blood Pressure (mmHg)': re.compile(r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:|measured)?\s*([\d/]+)(?:\s*mm\s*Hg)?', re.IGNORECASE),
                'Torsades de Pointes?': re.compile(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', re.IGNORECASE),
                'Medical History': re.compile(r'(?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
                'Medication History': re.compile(r'(?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
                'Course of Treatment': re.compile(r'(?:treatment|therapy|management|intervention|administered|given)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            }

            # Extract information using patterns
            for field, pattern in patterns.items():
                match = pattern.search(text)
                if match:
                    if field == 'Torsades de Pointes?':
                        report[field] = 'Yes'
                    else:
                        groups = [g for g in match.groups() if g]
                        report[field] = groups[0] if groups else ''

            return report

        except Exception as e:
            logger.error(f"Error analyzing paper: {str(e)}")
            return {}

    def _extract_with_context(self, text: str, pattern: re.Pattern) -> str:
        """Extract value and context using pattern."""
        if not text:
            return ""
            
        match = pattern.search(text)
        if not match:
            return ""
            
        # Get the matched text and its context
        matched_text = match.group(0)
        context = self.get_smart_context(text, match.start(), match.end())
        
        # Return the first group if available, otherwise the whole match
        groups = [g for g in match.groups() if g is not None]
        value = groups[0] if groups else matched_text
        
        return value

    def _extract_numeric(self, text: str) -> Optional[float]:
        """Extract numeric value from text."""
        if not text:
            return None
            
        try:
            # Find first number in text
            match = re.search(r'(\d+\.?\d*)', text)
            if match:
                return float(match.group(1))
        except Exception as e:
            logger.error(f"Error extracting numeric value: {e}")
            
        return None

    def _extract_unit(self, text: str) -> str:
        """Extract unit from dosing text."""
        if not text:
            return ""
            
        try:
            # Look for common units after numbers
            match = re.search(r'\d+\.?\d*\s*([a-zA-Z]+)', text)
            if match:
                return match.group(1)
        except Exception as e:
            logger.error(f"Error extracting unit: {e}")
            
        return ""

    def _extract_frequency(self, text: str) -> str:
        """Extract frequency from dosing text."""
        if not text:
            return ""
            
        try:
            # Look for common frequency patterns
            freq_patterns = [
                r'(once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid)',
                r'(per|a|each|every|q|q\.)\s*(?:day|daily|d|qd|od|q\.?d\.?|24\s*h(?:ours?)?|diem)'
            ]
            
            for pattern in freq_patterns:
                match = re.search(pattern, text, re.IGNORECASE)
                if match:
                    return match.group(0)
        except Exception as e:
            logger.error(f"Error extracting frequency: {e}")
            
        return ""

    def _clean_sex(self, text: str) -> str:
        """Clean and standardize sex value."""
        if not text:
            return ""
            
        text = text.lower()
        if any(term in text for term in ['male', 'man', 'm/', 'boy']):
            return "Male"
        elif any(term in text for term in ['female', 'woman', 'f/', 'girl']):
            return "Female"
            
        return ""

    def extract_value(self, text: str, pattern: re.Pattern) -> Optional[str]:
        """Extract value from text using a compiled regex pattern."""
        if not text:
            return None
            
        match = pattern.search(text)
        if match:
            # Return the first non-None group
            return next((g for g in match.groups() if g is not None), None)
        return None

    def extract_numeric(self, value: str) -> Optional[float]:
        """Extract numeric value from string."""
        if not value or not isinstance(value, str):
            return None
        match = re.search(r'(\d+\.?\d*)', value)
        return float(match.group(1)) if match else None

    def _extract_first_match(self, pattern, text):
        """Helper method to extract first regex match from text."""
        if not text:
            return None
            
        try:
            match = pattern.search(text)
            if match:
                # Return the first non-None group
                groups = [g for g in match.groups() if g is not None]
                return groups[0] if groups else match.group(0)
        except Exception as e:
            logger.error(f"Error in regex matching: {e}")
            
        return None

def convert_pdf_to_text(pdf_path: str) -> Optional[str]:
    """Convert a PDF file to text and save it"""
    try:
        # Create text file path
        text_path = pdf_path.rsplit('.', 1)[0] + '.txt'
        
        # Skip if text file already exists
        if os.path.exists(text_path):
            logger.info(f"Text file already exists: {text_path}")
            with open(text_path, 'r', encoding='utf-8') as f:
                return f.read()
            
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
                return text.strip()
                    
            except Exception as e:
                logger.error(f"Error reading PDF {pdf_path}: {e}")
                return None
                
    except Exception as e:
        logger.error(f"Error converting PDF to text: {e}")
        return None
