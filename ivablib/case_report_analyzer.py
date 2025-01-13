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
    text = ' '.join(text.split())
    text = text.replace('|', 'I').replace('—', '-').replace('–', '-')
    text = text.replace(''', "'").replace('"', '"').replace('"', '"')
    text = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', text)
    text = re.sub(r'(\d)\s+([.,])\s*(\d)', r'\1\2\3', text)
    return text

class CaseReportAnalyzer:
    """Analyzes medical case reports for relevant information."""

    def __init__(self):
        """Initialize the analyzer with regex patterns."""
        self.patterns = {
            'age': re.compile(r'(?:(?:aged?|age[d\s:]*|a|was)\s*)?(\d+)[\s-]*(?:year|yr|y|yo|years?)[s\s-]*(?:old|of\s+age)?|(?:age[d\s:]*|aged\s*)(\d+)', re.IGNORECASE),
            'sex': re.compile(r'\b(?:male|female|man|woman|boy|girl|[MF]\s*/\s*(?:\d+|[MF])|gender[\s:]+(?:male|female)|(?:he|she|his|her)\s+was)\b', re.IGNORECASE),
            'oral_dose': re.compile(r'(?:oral|po|by\s+mouth|orally|per\s+os)\s*(?:dose|dosage|doses|administration|administered|given|taken|prescribed|started|initiated|received|treated\s+with)?\s*(?:with|at|of|a|total|daily|once|twice|thrice|[1-4]x|q\.?d|bid|tid|qid)?\s*(?:dose\s+of\s+)?(\d+\.?\d*)\s*(?:mg|milligrams?|g|grams?|mcg|micrograms?|µg)', re.IGNORECASE),
            'qtc': re.compile(r'(?:QTc[FBR]?|corrected\s+QT|QT\s*corrected|QT[c]?\s*interval\s*(?:corrected)?|corrected\s*QT\s*interval)[\s:]*(?:interval|duration|measurement|prolongation|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
            'qt': re.compile(r'\b(?:QT|uncorrected\s+QT|QT\s*interval|interval\s*QT|baseline\s*QT)[\s:]*(?:interval|duration|measurement|value)?\s*(?:of|was|is|=|:|measured|found|documented|recorded|increased\s+to)?\s*(\d+\.?\d*)\s*(?:ms(?:ec)?|milliseconds?|s(?:ec)?|seconds?)?', re.IGNORECASE),
            'heart_rate': re.compile(r'(?:heart\s+rate|HR|pulse|ventricular\s+rate|heart\s+rhythm|cardiac\s+rate)[\s:]*(?:of|was|is|=|:|decreased\s+to|increased\s+to)?\s*(\d+)(?:\s*(?:beats?\s*per\s*min(?:ute)?|bpm|min-1|/min))?', re.IGNORECASE),
            'blood_pressure': re.compile(r'(?:blood\s+pressure|BP|arterial\s+pressure)[\s:]*(?:of|was|is|=|:|measured)?\s*([\d/]+)(?:\s*mm\s*Hg)?', re.IGNORECASE),
            'tdp': re.compile(r'\b(?:torsade[s]?\s*(?:de)?\s*pointes?|TdP|torsades|polymorphic\s+[vV]entricular\s+[tT]achycardia|PVT\s+(?:with|showing)\s+[tT]dP)\b', re.IGNORECASE),
            'medical_history': re.compile(r'(?:medical|clinical|past|previous|documented|known|significant)\s*(?:history|condition|diagnosis|comorbidities)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            'medication_history': re.compile(r'(?:medication|drug|prescription|current\s+medications?|concomitant\s+medications?)\s*(?:history|list|profile|regime)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE),
            'treatment_course': re.compile(r'(?:treatment|therapy|management|intervention|administered|given)[:\.]\s*([^\.]+?)(?:\.|\n|$)', re.IGNORECASE)
        }

    def analyze_papers(self, papers: List[Dict[str, Any]], columns: List[str]) -> pd.DataFrame:
        """Analyze multiple papers and return results as DataFrame."""
        results = []
        for paper in papers:
            try:
                result = self.analyze_paper(paper)
                if not result:
                    continue

                formatted_result = {
                    'Title': result['Title'],
                    'PMID': result['PMID'],
                    'DOI': result['DOI'],
                    'Abstract': result['Abstract'],
                    'Age': result['age_value'] or '',
                    'Sex': result['sex_value'] or '',
                    'Oral Dose (mg)': result['oral_dose_value'] or '',
                    'QTc (ms)': result['qtc_value'] or '',
                    'QT (ms)': result['qt_value'] or '',
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

        df = pd.DataFrame(results, columns=columns)
        df = df.replace({None: '', 'None': '', 'nan': ''})
        df = df.fillna('')
        return df

    def analyze_paper(self, paper: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Analyze a single paper for case report information."""
        try:
            text = paper.get('FullText', '') + ' ' + paper.get('Abstract', '')
            if not text.strip():
                return {}

            report = {
                'Title': paper.get('Title', ''),
                'PMID': paper.get('PMID', ''),
                'DOI': paper.get('DOI', ''),
                'Abstract': paper.get('Abstract', ''),
                'age_value': self._extract_numeric(text, self.patterns['age']),
                'sex_value': self._clean_sex(text, self.patterns['sex']),
                'oral_dose_value': self._extract_numeric(text, self.patterns['oral_dose']),
                'qtc_value': self._extract_numeric(text, self.patterns['qtc']),
                'qt_value': self._extract_numeric(text, self.patterns['qt']),
                'heart_rate_value': self._extract_numeric(text, self.patterns['heart_rate']),
                'blood_pressure_value': self._extract_first_match(text, self.patterns['blood_pressure']),
                'tdp_present': bool(self.patterns['tdp'].search(text)),
                'medical_history': self._extract_first_match(text, self.patterns['medical_history']),
                'medication_history': self._extract_first_match(text, self.patterns['medication_history']),
                'treatment_course': self._extract_first_match(text, self.patterns['treatment_course'])
            }
            return report
        except Exception as e:
            logger.error(f"Error analyzing paper: {e}")
            return None

    def _extract_numeric(self, text: str, pattern: re.Pattern) -> Optional[float]:
        """Extract numeric value from text using pattern."""
        if not text:
            return None
        match = pattern.search(text)
        if match:
            groups = [g for g in match.groups() if g is not None]
            if groups:
                try:
                    return float(groups[0])
                except (ValueError, TypeError):
                    pass
        return None

    def _clean_sex(self, text: str, pattern: re.Pattern) -> str:
        """Clean and standardize sex value."""
        if not text:
            return ""
        match = pattern.search(text)
        if match:
            value = match.group(0).lower()
            if any(term in value for term in ['male', 'man', 'm/', 'boy']):
                return "Male"
            elif any(term in value for term in ['female', 'woman', 'f/', 'girl']):
                return "Female"
        return ""

    def _extract_first_match(self, text: str, pattern: re.Pattern) -> str:
        """Extract first match from text using pattern."""
        if not text:
            return ""
        match = pattern.search(text)
        if match:
            groups = [g for g in match.groups() if g is not None]
            return groups[0] if groups else match.group(0)
        return ""
