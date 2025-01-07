"""Natural Language Processing for paper analysis."""
import re
from typing import Dict, List, Tuple, Union

def extract_qt_hr_values(text: str) -> Dict[str, List[Union[int, Tuple[int, str, str]]]]:
    """Extract QT, QTc, and HR values from text with context."""
    values = {
        'qt_values': [],
        'qtc_values': [],
        'hr_values': [],
        'bp_values': []
    }
    
    # QT interval patterns
    qt_patterns = [
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:of|was|is|=|:|measured|prolonged|increased|decreased)\s*(?:to|at|as)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:duration|measurement|value)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:baseline\s+)?QT\s+interval\s+(?:range|varied|ranging|from)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QT\s*(?:interval)?\s*(?:prolongation|increase|decrease)\s*(?:to|of|by)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})'
    ]
    
    # QTc patterns
    qtc_patterns = [
        r'(?P<pre>[^.]{0,100})QTc\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:of|was|is|=|:|measured)\s*(?:to|at|as)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})QTc\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:interval|duration|measurement|value)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:corrected|heart[\s-]rate[\s-]corrected)\s+QT\s*(?:\((?P<formula>Bazett|Fridericia|Framingham|Hodges)\))?\s*(?:was|is|=|:|of)?\s*(?P<value>\d{2,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:ms|msec|milliseconds)(?P<post>[^.]{0,100})'
    ]
    
    # Heart rate patterns
    hr_patterns = [
        r'(?P<pre>[^.]{0,100})(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:of|was|is|=|:|measured)\s*(?:to|at|as)?\s*(?P<value>\d{1,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:increased|decreased|changed)\s*(?:to|at)?\s*(?P<value>\d{1,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})',
        r'(?P<pre>[^.]{0,100})(?:mean|average)\s+(?:heart\s+rate|HR|pulse(?:\s+rate)?)\s*(?:was|is|=|:|of)?\s*(?P<value>\d{1,3}(?:\.\d)?)\s*(?:±\s*\d+(?:\.\d)?)?\s*(?:bpm|beats\s*(?:per|\/)?\s*min(?:ute)?)(?P<post>[^.]{0,100})'
    ]
    
    def extract_with_context(patterns, text, value_range, value_type):
        results = []
        for pattern in patterns:
            matches = re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE | re.DOTALL)
            for match in matches:
                try:
                    value = float(match.group('value'))
                    min_val, max_val = value_range
                    
                    if min_val <= value <= max_val:
                        pre_context = match.group('pre').strip()
                        post_context = match.group('post').strip()
                        context = f"{pre_context} [...] {post_context}"
                        context = re.sub(r'\s+', ' ', context).strip()
                        
                        if value_type == 'qtc' and 'formula' in match.groupdict():
                            formula = match.group('formula') or 'unspecified'
                            results.append((int(value), formula, context))
                        else:
                            method = 'measured' if any(word in pre_context.lower() 
                                for word in ['measured', 'recorded', 'observed', 'mean', 'average']) else 'reported'
                            results.append((int(value), method, context))
                except (ValueError, IndexError):
                    continue
        return results
    
    # Extract values with appropriate ranges
    values['qt_values'] = extract_with_context(qt_patterns, text, (200, 600), 'qt')
    values['qtc_values'] = extract_with_context(qtc_patterns, text, (200, 600), 'qtc')
    values['hr_values'] = extract_with_context(hr_patterns, text, (30, 200), 'hr')
    
    return values

def extract_case_details(text: str) -> Dict[str, str]:
    """Extract detailed case information using NLP patterns."""
    details = {
        'age': None,
        'sex': None,
        'medical_history': [],
        'medications': [],
        'outcome': None
    }
    
    # Age patterns
    age_patterns = [
        r'(\d+)[-\s]year[-\s]old',
        r'age[d]?\s+(\d+)',
        r'(\d+)\s*(?:yo|y\.o\.|years\s+old)'
    ]
    
    for pattern in age_patterns:
        match = re.search(pattern, text.lower())
        if match:
            try:
                age = int(match.group(1))
                if 0 <= age <= 120:
                    details['age'] = age
                    break
            except ValueError:
                continue
    
    # Sex patterns
    sex_match = re.search(r'\b(male|female|man|woman|boy|girl)\b', text.lower())
    if sex_match:
        sex = sex_match.group(1)
        details['sex'] = 'Male' if sex in ['male', 'man', 'boy'] else 'Female'
    
    # Medical history patterns
    history_patterns = [
        r'(?:medical|past) history (?:of|significant for|notable for) (.*?)\.',
        r'(?:diagnosed with|known) (.*?) (?:prior to|before)',
        r'(?:history of) (.*?) (?:was|is|noted)'
    ]
    
    for pattern in history_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            history = match.group(1).strip()
            if history and len(history) < 100:  # Reasonable length check
                details['medical_history'].append(history)
    
    return details
