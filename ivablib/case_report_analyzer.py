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
            # Enhanced demographic patterns with more variations
            'age': re.compile(r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
            
            'sex': re.compile(r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b', re.IGNORECASE),
            
            # Enhanced dosing patterns with medical abbreviations
            'oral_dose': re.compile(r"""
                (?:
                    (?:oral|po|by\s+mouth|orally|per\s+os)\s+
                    (?:dose|dosage|doses|administration|administered|given|taken|prescribed|started|initiated|received|treated\s+with)|
                    (?:dose[d\s]*|given|administered|took|ingested|consumed|prescribed|started|initiated|received|treated\s+with)\s+
                    (?:orally|po|by\s+mouth|per\s+os)
                )
                \s*
                (?:with|at|of|a|total|daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid)?\s*
                (?:dose\s+of\s+)?
                (\d+\.?\d*)\s*
                (?:mg|milligrams?|g|grams?|mcg|micrograms?|µg)
                (?:\s*(?:per|a|each|every|q|q\.)\s*(?:day|daily|d|qd|od|q\.?d\.?|24\s*h(?:ours?)?|diem))?
            """, re.VERBOSE | re.IGNORECASE),
            
            # Enhanced ECG patterns with more variations
            'qtc': re.compile(r"""
                (?:
                    QTc[FBR]?[\s-]*(?:interval|duration|measurement|value)?|
                    corrected\s+QT[\s-]*(?:interval|duration|measurement|value)?|
                    QT[\s-]*(?:interval|duration|measurement|value)?\s*corrected|
                    QT[c]?[\s-]*interval[\s-]*(?:corrected|correction)|
                    corrected[\s-]*QT[\s-]*interval
                )
                [\s:=-]*
                (?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to|prolonged\s+to|approximately|~|≈|≅|about)?
                \s*
                (\d+(?:\.\d+)?)\s*
                (?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?
                (?:\s*\([^)]*\))?  # Capture optional parenthetical notes
            """, re.VERBOSE | re.IGNORECASE),

            'qt': re.compile(r"""
                (?:
                    (?<!QTc)(?<!QTcF)(?<!QTcB)(?<!QTcR)  # Negative lookbehind to avoid matching QTc
                    (?:
                        QT[\s-]*interval|
                        QT[\s-]*duration|
                        QT[\s-]*measurement|
                        QT[\s-]*prolongation|
                        QT[\s-]*value|
                        QT[\s-]*=|
                        QT[\s-]*:|
                        QT[\s-]*was|
                        QT[\s-]*is|
                        QT[\s-]*of|
                        \bQT\b
                    )
                )
                \s*
                (?:interval|duration|measurement|value)?
                [\s:=-]*
                (?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to|prolonged\s+to|approximately|~|≈|≅|about)?
                \s*
                (\d+(?:\.\d+)?)\s*
                (?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?
                (?:\s*\([^)]*\))?  # Capture optional parenthetical notes
            """, re.VERBOSE | re.IGNORECASE),
            
            'heart_rate': re.compile(r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:|decreased\s+to|increased\s+to)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?', re.IGNORECASE),
            
            'blood_pressure': re.compile(r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:|measured)?\s*([\d/]+)(?:\s*mm\s*Hg)?', re.IGNORECASE),
            
            # Enhanced arrhythmia patterns
            'tdp': re.compile(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', re.IGNORECASE),
            
            # Enhanced history patterns with better context capture
            'medical_history': re.compile(r"""
                (?:
                    (?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)|
                    (?:history\s+of|diagnosed\s+with|presented\s+with|admitted\s+with|suffering\s+from)|
                    (?:comorbid\s+conditions?|underlying\s+conditions?|pre-existing\s+conditions?)
                )
                [:\.;\s]+
                ([^\.]+?)
                (?:\.(?:\s+[A-Z]|\n|$)|$)
            """, re.VERBOSE | re.IGNORECASE),
            
            'medication_history': re.compile(r"""
                (?:
                    (?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)|
                    (?:medications?\s+included|taking|receiving|treated\s+with|on\s+the\s+following\s+medications?)|
                    (?:home\s+medications?|outpatient\s+medications?|chronic\s+medications?)
                )
                [:\.;\s]+
                ([^\.]+?)
                (?:\.(?:\s+[A-Z]|\n|$)|$)
            """, re.VERBOSE | re.IGNORECASE),
            
            'treatment_course': re.compile(r"""
                (?:
                    (?:treatment|therapy|management|intervention|administered|given)|
                    (?:treated|managed|received|started\s+on|initiated\s+on)|
                    (?:course\s+of\s+treatment|therapeutic\s+approach|clinical\s+course)
                )
                [:\.;\s]+
                ([^\.]+?)
                (?:\.(?:\s+[A-Z]|\n|$)|$)
            """, re.VERBOSE | re.IGNORECASE),
            
            # Enhanced treatment patterns
            'drug_qt': re.compile(r"""
                (?:
                    (?:QT\s+prolongation|QT\s+interval\s+prolongation|QTc\s+prolongation)\s+
                    (?:associated\s+with|due\s+to|caused\s+by|induced\s+by|related\s+to)\s+
                    (?:the\s+use\s+of|the\s+administration\s+of|the\s+prescription\s+of)\s+
                    ([\w\s,]+)
                )
            """, re.VERBOSE | re.IGNORECASE),
            
            'brady_qt': re.compile(r"""
                (?:
                    (?:bradycardia|bradyarrhythmia)\s+
                    (?:associated\s+with|due\s+to|caused\s+by|induced\s+by|related\s+to)\s+
                    (?:QT\s+prolongation|QT\s+interval\s+prolongation|QTc\s+prolongation)
                )
            """, re.VERBOSE | re.IGNORECASE),
            
            'figure_qt': re.compile(r"""
                (?:
                    (?:figure|fig\.?)\s+
                    [\d\.]+\s*
                    (?:shows|demonstrates|depicts|displays)\s+
                    (?:QT\s+prolongation|QT\s+interval\s+prolongation|QTc\s+prolongation)
                )
            """, re.VERBOSE | re.IGNORECASE),
            
            'table_qt': re.compile(r"""
                (?:
                    (?:table|tab\.?)\s+
                    [\d\.]+\s*
                    (?:lists|summarizes|presents)\s+
                    (?:QT\s+prolongation|QT\s+interval\s+prolongation|QTc\s+prolongation)
                )
            """, re.VERBOSE | re.IGNORECASE),
            
            'supp_qt': re.compile(r"""
                (?:
                    (?:supplementary|supplement)\s+
                    (?:material|data|information)\s+
                    (?:includes|contains|provides)\s+
                    (?:QT\s+prolongation|QT\s+interval\s+prolongation|QTc\s+prolongation)
                )
            """, re.VERBOSE | re.IGNORECASE),
        }

    def get_smart_context(self, text: str, match_start: int, match_end: int, min_chars: int = 200) -> str:
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

    def analyze(self, text: str) -> Dict[str, Any]:
        """
        Analyze text to extract relevant information.
        
        Args:
            text: The text to analyze
            
        Returns:
            Dict containing extracted information
        """
        # Clean text
        text = self._clean_text(text)
        
        # Initialize results
        results = {
            'age': '',
            'sex': '',
            'oral_dose_value': '',
            'oral_dose_unit': '',
            'oral_dose_freq': '',
            'qt_value': '',
            'qtc_value': '',
            'heart_rate_value': '',
            'blood_pressure_value': '',
            'tdp_present': False,
            'tdp_context': '',
            'medical_history': '',
            'medication_history': '',
            'treatment_course': ''
        }
        
        # Extract each field with context
        for field, pattern in self.patterns.items():
            matches = list(pattern.finditer(text))
            if matches:
                if field == 'tdp':
                    results['tdp_present'] = True
                    # Get context for the first match
                    match = matches[0]
                    results['tdp_context'] = self.get_smart_context(text, match.start(), match.end())
                elif field == 'oral_dose':
                    # Use the first match for dosing info
                    match = matches[0]
                    dose_text = match.group(0)
                    value_match = re.search(r'(\d+(?:\.\d+)?)', dose_text)
                    unit_match = re.search(r'(mg|milligrams?|g|grams?|mcg|micrograms?|µg)', dose_text, re.IGNORECASE)
                    freq_match = re.search(r'(daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid|per\s+day|/day)', dose_text, re.IGNORECASE)
                    
                    results['oral_dose_value'] = value_match.group(1) if value_match else ''
                    results['oral_dose_unit'] = unit_match.group(1) if unit_match else ''
                    results['oral_dose_freq'] = freq_match.group(1) if freq_match else ''
                elif field in ['qt', 'qtc']:
                    # For QT measurements, try to find the most relevant value
                    # Sort matches by the numeric value to find the most concerning one
                    values = []
                    for match in matches:
                        try:
                            # Look for the numeric value in all groups
                            numeric_groups = [g for g in match.groups() if g and re.match(r'^\d+(?:\.\d+)?$', g)]
                            if numeric_groups:
                                value = float(numeric_groups[0])
                                if 100 <= value <= 700:  # Reasonable QT range in ms
                                    context = self.get_smart_context(text, match.start(), match.end())
                                    values.append((value, context))
                        except (ValueError, IndexError):
                            continue
                    
                    if values:
                        # For QTc, prefer larger values as they're more concerning
                        # For QT, also prefer larger values
                        value, context = max(values, key=lambda x: x[0])
                        results[f'{field}_value'] = str(value)
                elif field in ['heart_rate', 'blood_pressure']:
                    # For numeric fields, get the first numeric value
                    match = matches[0]
                    for group in match.groups():
                        if group and re.match(r'^\d+(?:/\d+)?$', group):  # Allow for BP format (e.g., 120/80)
                            results[f'{field}_value'] = group.strip()
                            break
                elif field in ['medical_history', 'medication_history', 'treatment_course']:
                    # For history fields, get the content after the pattern
                    match = matches[0]
                    if match.groups():
                        content = next((g for g in match.groups() if g), '').strip()
                        results[field] = content
                elif field == 'sex':
                    # Clean up sex field
                    match = matches[0]
                    sex_text = match.group(0).lower()
                    if any(term in sex_text for term in ['male', 'man', 'boy', 'm/']):
                        results['sex'] = 'Male'
                    elif any(term in sex_text for term in ['female', 'woman', 'girl', 'f/']):
                        results['sex'] = 'Female'
                elif field == 'age':
                    # Get the first numeric age value
                    match = matches[0]
                    for group in match.groups():
                        if group and group.isdigit():
                            results['age'] = group
                            break
    
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
        
    def analyze_paper(self, paper: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Analyze a single paper for case report information."""
        try:
            # Get text from both full text and abstract
            text = ' '.join(filter(None, [
                paper.get('FullText', ''),
                paper.get('Abstract', '')
            ]))
            
            if not text.strip():
                logger.warning(f"No text content found for paper {paper.get('PMID', 'Unknown')}")
                return {}

            # Clean text
            text = self._clean_text(text)
            
            # Extract values with context
            result = {
                'title': paper.get('Title', ''),
                'pmid': paper.get('PMID', ''),
                'doi': paper.get('DOI', '')
            }

            # Extract demographic info
            age_match = self._extract_with_context(text, self.patterns['age'])
            result['age'] = age_match[0] if age_match else None
            result['age_context'] = age_match[1] if age_match else None

            sex_match = self._extract_with_context(text, self.patterns['sex'])
            result['sex'] = self._clean_sex(sex_match[0]) if sex_match else None
            result['sex_context'] = sex_match[1] if sex_match else None

            # Extract dosing info
            dose_match = self._extract_with_context(text, self.patterns['oral_dose'])
            if dose_match:
                result['oral_dose_value'] = self._extract_numeric(dose_match[0])
                result['oral_dose_unit'] = self._extract_unit(dose_match[0])
                result['oral_dose_freq'] = self._extract_frequency(dose_match[0])
                result['oral_dose_context'] = dose_match[1]

            # Extract ECG measurements
            qt_match = self._extract_with_context(text, self.patterns['qt'])
            if qt_match:
                result['qt_value'] = self._extract_numeric(qt_match[0])
                result['qt_context'] = qt_match[1]

            qtc_match = self._extract_with_context(text, self.patterns['qtc'])
            if qtc_match:
                result['qtc_value'] = self._extract_numeric(qtc_match[0])
                result['qtc_context'] = qtc_match[1]

            # Extract vital signs
            hr_match = self._extract_with_context(text, self.patterns['heart_rate'])
            if hr_match:
                result['heart_rate_value'] = self._extract_numeric(hr_match[0])
                result['heart_rate_context'] = hr_match[1]

            bp_match = self._extract_with_context(text, self.patterns['blood_pressure'])
            if bp_match:
                result['blood_pressure_value'] = bp_match[0]
                result['blood_pressure_context'] = bp_match[1]

            # Check for TdP
            tdp_match = self._extract_with_context(text, self.patterns['tdp'])
            result['tdp_present'] = bool(tdp_match)
            result['tdp_context'] = tdp_match[1] if tdp_match else None

            # Extract history
            med_history = self._extract_with_context(text, self.patterns['medical_history'])
            result['medical_history'] = med_history[0] if med_history else None
            result['medical_history_context'] = med_history[1] if med_history else None

            med_list = self._extract_with_context(text, self.patterns['medication_history'])
            result['medication_history'] = med_list[0] if med_list else None
            result['medication_history_context'] = med_list[1] if med_list else None

            treatment = self._extract_with_context(text, self.patterns['treatment_course'])
            result['treatment_course'] = treatment[0] if treatment else None
            result['treatment_course_context'] = treatment[1] if treatment else None

            return result

        except Exception as e:
            logger.error(f"Error analyzing paper {paper.get('PMID', 'Unknown')}: {str(e)}")
            return {}

    def _extract_with_context(self, text: str, pattern: re.Pattern) -> Optional[Tuple[str, str]]:
        """Extract value and context using pattern."""
        if not text:
            return None
            
        match = pattern.search(text)
        if not match:
            return None
            
        value = next((g for g in match.groups() if g is not None), None)
        if not value:
            return None
            
        context = self.get_smart_context(text, match.start(), match.end())
        return value, context

    def _extract_numeric(self, text: str) -> Optional[float]:
        """Extract numeric value from text."""
        if not text:
            return None
            
        match = re.search(r'(\d+\.?\d*)', text)
        return float(match.group(1)) if match else None

    def _extract_unit(self, text: str) -> Optional[str]:
        """Extract unit from dosing text."""
        if not text:
            return None
            
        unit_pattern = re.compile(r'(mg|milligrams?|g|grams?|mcg|micrograms?|µg)', re.IGNORECASE)
        match = unit_pattern.search(text)
        return match.group(1) if match else None

    def _extract_frequency(self, text: str) -> Optional[str]:
        """Extract frequency from dosing text."""
        if not text:
            return None
            
        freq_pattern = re.compile(r'(daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid|per\s+day|/day)', re.IGNORECASE)
        match = freq_pattern.search(text)
        return match.group(1) if match else None

    def _clean_sex(self, text: str) -> Optional[str]:
        """Clean and standardize sex value."""
        if not text:
            return None
            
        text = text.lower()
        if any(m in text for m in ['male', 'man', 'm/', 'boy']):
            return 'Male'
        elif any(f in text for f in ['female', 'woman', 'f/', 'girl']):
            return 'Female'
        return None

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

    def _clean_text(self, text: str) -> str:
        """Clean and normalize text for analysis."""
        if not text:
            return ""
            
        # Preserve case for important medical terms
        text = re.sub(r'\bQT[cC]?\b', lambda m: m.group(0).upper(), text)
        text = re.sub(r'\b[Tt]d[Pp]\b', 'TdP', text)
        text = re.sub(r'\b[Ee][Cc][Gg]\b', 'ECG', text)
        text = re.sub(r'\b[Ee][Kk][Gg]\b', 'ECG', text)
        
        # Replace common medical abbreviations
        replacements = {
            'yo ': ' year old ',
            'y/o': ' year old ',
            'yr ': ' year ',
            'yrs': ' years ',
            'yo,': ' year old,',
            'y/o,': ' year old,',
            'yr,': ' year,',
            'yrs,': ' years,',
            'hr ': ' heart rate ',
            'bp ': ' blood pressure ',
            'hx': ' history ',
            'dx': ' diagnosis ',
            'tx': ' treatment ',
            'pmh': ' past medical history ',
            'psh': ' past surgical history ',
            'meds': ' medications ',
            'sig': ' significant ',
            'w/': ' with ',
            'w/o': ' without ',
            's/p': ' status post ',
            'c/o': ' complains of ',
            'h/o': ' history of ',
            'p/w': ' presented with ',
        }
        
        # Apply replacements
        for old, new in replacements.items():
            text = re.sub(r'\b' + re.escape(old) + r'\b', new, text, flags=re.IGNORECASE)
            
        # Remove excessive whitespace
        text = re.sub(r'\s+', ' ', text)
        
        # Normalize sentence endings
        text = re.sub(r'([.!?])\s*([A-Z])', r'\1 \2', text)
        
        return text.strip()

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
