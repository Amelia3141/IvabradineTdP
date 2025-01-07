import unittest
from ivablib.case_report_analyzer import CaseReportAnalyzer

class TestCaseReportAnalyzer(unittest.TestCase):
    def setUp(self):
        self.analyzer = CaseReportAnalyzer()
        
    def test_extract_age(self):
        """Test age extraction from text"""
        test_cases = [
            ("A 73-year-old woman presented with...", 73),
            ("The patient was an 85-year-old man", 85),
            ("Age: 80 years", 80),
            ("Patient aged 83", 83),
        ]
        
        for text, expected_age in test_cases:
            case = {"title": "", "full_text": text}
            self.analyzer._extract_info(case)
            self.assertEqual(case.get("age"), expected_age)
            
    def test_extract_sex(self):
        """Test sex extraction from text"""
        test_cases = [
            ("A 73-year-old woman presented with...", "woman"),
            ("The patient was an 85-year-old man", "man"),
            ("Female patient aged 83", "woman"),
            ("Male, age 80 years", "man"),
        ]
        
        for text, expected_sex in test_cases:
            case = {"title": "", "full_text": text}
            self.analyzer._extract_info(case)
            self.assertEqual(case.get("sex"), expected_sex)

if __name__ == '__main__':
    unittest.main()
