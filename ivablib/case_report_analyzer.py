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

class CaseReportAnalyzer:
    """Analyzes medical case reports for relevant information."""
    
    def __init__(self):
        """Initialize the analyzer with regex patterns."""
        # Compile regex patterns for better performance
        self.patterns = {
            # Enhanced demographic patterns
            'age': re.compile(r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
            
            'sex': re.compile(r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b', re.IGNORECASE),
            
            # Enhanced dosing patterns
            'oral_dose': re.compile(r'(?:oral\s+dose|dose[d\s]*orally|given|administered|took|ingested|consumed|prescribed)\s*(?:with|at|of|a|total)?\s*(\d+\.?\d*)\s*(?:mg|milligrams?|g|grams?|mcg|micrograms?)', re.IGNORECASE),
            
            # Enhanced concentration patterns
            'theoretical_max': re.compile(r'(?:maximum|max|peak|highest|cmax|c-?max)\s*(?:concentration|conc|level|exposure)\s*(?:of|was|is|reached|measured|observed|calculated)?\s*(?:approximately|about|~|≈)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar|nmol/[lL]|μmol/[lL]|ng/m[lL])', re.IGNORECASE),
            
            'bioavailability': re.compile(r'(?:bioavailability|F|absorption|systemic\s+availability)\s*(?:of|was|is|approximately|about|~)?\s*(?:approximately|about|~|≈)?\s*(\d+\.?\d*)\s*%?', re.IGNORECASE),
            
            'herg_ic50': re.compile(r'(?:hERG\s+IC50|IC50\s+(?:value\s+)?(?:for|of)\s+hERG|half-?maximal\s+(?:inhibitory)?\s+concentration|IC50)\s*(?:of|was|is|=|:)?\s*(?:approximately|about|~|≈)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar|nmol/[lL]|μmol/[lL])', re.IGNORECASE),
            
            'plasma_concentration': re.compile(r'(?:plasma|blood|serum)\s*(?:concentration|level|exposure|Cmax|C-max)\s*(?:of|was|is|measured|found|detected|reached)?\s*(?:approximately|about|~|≈)?\s*(\d+\.?\d*)\s*(?:μM|uM|micromolar|ng/m[lL]|μg/m[lL]|mg/[lL]|g/[lL])', re.IGNORECASE),
            
            # Enhanced ECG patterns
            'qtc': re.compile(r"""
                (?:
                    QTc[FBR]?|
                    corrected\s+QT|
                    QT\s*corrected|
                    QT[c]?\s*interval\s*(?:corrected)?|
                    corrected\s*QT\s*interval
                )
                [\s:]*
                (?:interval|duration|measurement|prolongation|value)?
                \s*
                (?:of|was|is|=|:|measured|found|documented|recorded)?
                \s*
                (\d+\.?\d*)
                \s*
                (?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?
            """, re.VERBOSE | re.IGNORECASE),

            'qt': re.compile(r"""
                \b(?:
                    QT|
                    uncorrected\s+QT|
                    QT\s*interval|
                    interval\s*QT|
                    baseline\s*QT
                )
                [\s:]*
                (?:interval|duration|measurement|value)?
                \s*
                (?:of|was|is|=|:|measured|found|documented|recorded)?
                \s*
                (\d+\.?\d*)
                \s*
                (?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?
            """, re.VERBOSE | re.IGNORECASE),
            
            'heart_rate': re.compile(r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?', re.IGNORECASE),
            
            'blood_pressure': re.compile(r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:)?\s*([\d/]+)(?:\s*mm\s*Hg)?', re.IGNORECASE),
            
            # Enhanced arrhythmia patterns
            'tdp': re.compile(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', re.IGNORECASE),
            
            # Enhanced history patterns
            'medical_history': re.compile(r'(?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            'medication_history': re.compile(r'(?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            # Enhanced treatment patterns
            'treatment_course': re.compile(r'(?:treatment|therapy|management|intervention|administered|given)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            # Additional patterns for table/figure data
            'table_data': re.compile(r'Table\s*\d+[.:]\s*([^.]+)', re.IGNORECASE),
            
            'figure_data': re.compile(r'Fig(?:ure|\.)\s*\d+[.:]\s*([^.]+)', re.IGNORECASE),
            
            # Laboratory values
            'lab_values': re.compile(r'(?:laboratory|lab|biochemistry|chemistry)\s*(?:results|values|findings|data)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            
            # Figure/Table QT mentions
            'figure_qt': re.compile(r"""
                (?:Fig(?:ure|\.)?|Table)\s*\d+[A-Z]?\)?[:\s]+
                (?:[^.]*?)
                (?:
                    QT[c]?[\s:]*(\d+\.?\d*)(?:\s*(?:ms|msec|s|sec|seconds?))?|
                    (?:shows?|demonstrates?|illustrates?|reveals?)\s+
                    (?:[^.]*?)
                    QT[c]?[\s:]*(\d+\.?\d*)
                )
            """, re.VERBOSE | re.IGNORECASE),

            # Table cell QT values
            'table_qt': re.compile(r"""
                (?:QT[c]?|Interval)[\s:]*
                (?:\((?:ms|msec|s|sec|seconds?)\))?\s*
                [\|\{\[\t]*
                (\d+\.?\d*)
            """, re.VERBOSE | re.IGNORECASE),
            
            # Supplementary material QT values
            'supp_qt': re.compile(r"""
                (?:Supplementary|Suppl\.?|Online)\s+
                (?:Fig(?:ure|\.)?|Table|Data|Material)\s*\d*[:\s]+
                (?:[^.]*?)
                QT[c]?[\s:]*(\d+\.?\d*)
            """, re.VERBOSE | re.IGNORECASE),
            
            # Drug-QT relationships
            'drug_qt': re.compile(r"""
                (?P<drug>\w+)\s+
                (?:
                    (?:is|are|was|were|being|been)\s+
                    (?:
                        known\s+to|
                        reported\s+to|
                        associated\s+with|
                        linked\s+to
                    )\s+
                    (?:cause|produce|induce|prolong|affect)
                )?\s*
                (?:QT|QTc)\s*
                (?:
                    prolongation|
                    interval|
                    changes?
                )
            """, re.VERBOSE | re.IGNORECASE),
            
            # Bradycardia-QT relationship
            'brady_qt': re.compile(r"""
                (?:
                    bradycardia|
                    slow\s+heart\s+rate|
                    decreased\s+heart\s+rate|
                    heart\s+rate\s*(?:of)?\s*\d+
                )
                [^.]*?
                (?:QT|QTc|prolongation|interval)
            """, re.VERBOSE | re.IGNORECASE)
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
            return next((g for g in match.groups() if g is not None), None)
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
                text = paper.get('FullText', '') + ' ' + paper.get('Abstract', '')
                if not text.strip():
                    continue
                    
                # Extract all values with context
                result = {
                    'title': paper.get('Title', ''),
                    'age': self._extract_first_match(self.patterns['age'], text),
                    'sex': self._extract_first_match(self.patterns['sex'], text),
                    'oral_dose': self._extract_first_match(self.patterns['oral_dose'], text),
                    'theoretical_max': self._extract_first_match(self.patterns['theoretical_max'], text),
                    'bioavailability': self._extract_first_match(self.patterns['bioavailability'], text),
                    'herg_ic50': self._extract_first_match(self.patterns['herg_ic50'], text),
                    'plasma_concentration': self._extract_first_match(self.patterns['plasma_concentration'], text),
                    'qt': self._extract_first_match(self.patterns['qt'], text),
                    'qtc': self._extract_first_match(self.patterns['qtc'], text),
                    'heart_rate': self._extract_first_match(self.patterns['heart_rate'], text),
                    'tdp': 'Yes' if self.patterns['tdp'].search(text) else 'No',
                    'blood_pressure': self._extract_first_match(self.patterns['blood_pressure'], text),
                    'medical_history': self._extract_first_match(self.patterns['medical_history'], text),
                    'medication_history': self._extract_first_match(self.patterns['medication_history'], text),
                    'treatment_course': self._extract_first_match(self.patterns['treatment_course'], text)
                }
                
                # Calculate QTR and QTF if possible
                qt_val = self.extract_numeric(result['qt'])
                hr_val = self.extract_numeric(result['heart_rate'])
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
                    'Oral Dose (mg)': result['oral_dose'] or '',
                    'theoretical max concentration (μM)': result['theoretical_max'] or '',
                    '40% bioavailability': result['bioavailability'] or '',
                    'Theoretical HERG IC50 / Concentration μM': result['herg_ic50'] or '',
                    '40% Plasma concentration': result['plasma_concentration'] or '',
                    'Uncorrected QT (ms)': result['qt'] or '',
                    'QTc': result['qtc'] or '',
                    'QTR': result['qtr'] or '',
                    'QTF': result['qtf'] or '',
                    'Heart Rate (bpm)': result['heart_rate'] or '',
                    'Torsades de Pointes?': result['tdp'],
                    'Blood Pressure (mmHg)': result['blood_pressure'] or '',
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
        
    def extract_numeric(self, value: str) -> Optional[float]:
        """Extract numeric value from string."""
        if not value or not isinstance(value, str):
            return None
        match = re.search(r'(\d+\.?\d*)', value)
        return float(match.group(1)) if match else None

    def analyze_paper(self, paper: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Analyze a single paper for case report information."""
        try:
            # Get text from full text if available, otherwise use title and abstract
            text = paper.get('FullText', '')
            if not text:
                text = f"{paper.get('Title', '')} {paper.get('Abstract', '')}"
            
            # Skip if no meaningful text to analyze
            if not text.strip():
                return None
            
            # Perform comprehensive QT analysis
            qt_analysis = self.analyze_qt_comprehensive(text)
            
            # Extract QT/QTc values from comprehensive analysis
            qtc_value = None
            qt_value = None
            qtc_context = None
            qt_context = None
            qt_source = None
            
            # Function to extract numeric value from text
            def extract_numeric(text):
                match = re.search(r'(\d+\.?\d*)', text)
                return float(match.group(1)) if match else None
            
            # Check all sources for QT/QTc values in priority order
            sources = [
                ('Main Text', qt_analysis['numeric_values']),
                ('Figure', qt_analysis['figure_qt']),
                ('Table', qt_analysis['table_qt']),
                ('Supplementary', qt_analysis['supp_qt'])
            ]
            
            for source_name, values in sources:
                # For numeric values from main text
                if source_name == 'Main Text':
                    for pattern_name, pattern_values in values.items():
                        if pattern_values:
                            max_value_dict = max(pattern_values, key=lambda x: x['value'])
                            if 'qtc' in pattern_name and (qtc_value is None or max_value_dict['value'] > qtc_value):
                                qtc_value = max_value_dict['value']
                                qtc_context = max_value_dict['context']
                                qt_source = 'Main Text'
                            elif 'qt' in pattern_name and (qt_value is None or max_value_dict['value'] > qt_value):
                                qt_value = max_value_dict['value']
                                qt_context = max_value_dict['context']
                                qt_source = 'Main Text'
                else:
                    # For figure/table/supplementary values
                    for item in values:
                        value = extract_numeric(item['text'])
                        if value:
                            # Assume milliseconds if no unit specified
                            if 'sec' in item['text'].lower():
                                value *= 1000
                            
                            # Update QT/QTc based on context
                            if 'qtc' in item['text'].lower():
                                if qtc_value is None or value > qtc_value:
                                    qtc_value = value
                                    qtc_context = item['context']
                                    qt_source = source_name
                            else:
                                if qt_value is None or value > qt_value:
                                    qt_value = value
                                    qt_context = item['context']
                                    qt_source = source_name
            
            # Extract other information with context
            age_info = self.extract_value_with_context(text, self.patterns['age'])
            sex_info = self.extract_value_with_context(text, self.patterns['sex'])
            dose_info = self.extract_value_with_context(text, self.patterns['oral_dose'])
            max_conc_info = self.extract_value_with_context(text, self.patterns['theoretical_max'])
            bioavail_info = self.extract_value_with_context(text, self.patterns['bioavailability'])
            herg_info = self.extract_value_with_context(text, self.patterns['herg_ic50'])
            plasma_info = self.extract_value_with_context(text, self.patterns['plasma_concentration'])
            hr_info = self.extract_value_with_context(text, self.patterns['heart_rate'])
            bp_info = self.extract_value_with_context(text, self.patterns['blood_pressure'])
            
            result = {
                'title': paper.get('Title', ''),
                'year': paper.get('Year', ''),
                'abstract': paper.get('Abstract', ''),
                'full_text': paper.get('FullText', 'None'),
                'age': age_info[0] if age_info else None,
                'age_context': age_info[1] if age_info else None,
                'sex': sex_info[0] if sex_info else None,
                'sex_context': sex_info[1] if sex_info else None,
                'oral_dose': dose_info[0] if dose_info else None,
                'dose_context': dose_info[1] if dose_info else None,
                'theoretical_max': max_conc_info[0] if max_conc_info else None,
                'max_conc_context': max_conc_info[1] if max_conc_info else None,
                'bioavailability': bioavail_info[0] if bioavail_info else None,
                'bioavailability_context': bioavail_info[1] if bioavail_info else None,
                'herg_ic50': herg_info[0] if herg_info else None,
                'herg_context': herg_info[1] if herg_info else None,
                'plasma_concentration': plasma_info[0] if plasma_info else None,
                'plasma_context': plasma_info[1] if plasma_info else None,
                'qt': qt_value,
                'qt_context': qt_context,
                'qt_source': qt_source,
                'qtc': qtc_value,
                'qtc_context': qtc_context,
                'qtc_source': qt_source,
                'heart_rate': hr_info[0] if hr_info else None,
                'heart_rate_context': hr_info[1] if hr_info else None,
                'tdp': 'Yes' if qt_analysis['tdp_mentions'] else 'No',
                'tdp_context': qt_analysis['tdp_mentions'][0]['context'] if qt_analysis['tdp_mentions'] else None,
                'blood_pressure': bp_info[0] if bp_info else None,
                'blood_pressure_context': bp_info[1] if bp_info else None,
                'medical_history': self.extract_value(text, self.patterns['medical_history']),
                'medication_history': self.extract_value(text, self.patterns['medication_history']),
                'treatment_course': self.extract_value(text, self.patterns['treatment_course'])
            }
            
            # Add detailed QT analysis
            result['detailed_qt_analysis'] = qt_analysis
            
            # Only return if we found some meaningful information
            if any(v for v in result.values() if v not in (None, '', 'No')):
                return result
            return None
            
        except Exception as e:
            logger.error(f"Error analyzing paper: {str(e)}")
            return None

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

def analyze_papers(papers: List[Dict[str, Any]], drug_name: str) -> pd.DataFrame:
    """Analyze a list of papers and create a case report table."""
    analyzer = CaseReportAnalyzer()
    return analyzer.analyze_papers(papers, drug_name)

def clean_text(text: str) -> str:
    """Clean and normalize text extracted from PDF."""
    if not text:
        return ""
        
    # Remove excessive whitespace
    text = re.sub(r'\s+', ' ', text)
    
    # Fix common PDF extraction issues
    text = text.replace('- ', '')  # Remove hyphenation
    text = text.replace('•', '- ')  # Convert bullets to dashes
    
    return text.strip()

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
