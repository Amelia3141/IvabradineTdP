"""User interface handling."""
from typing import Optional, Dict, Any
import sys

def get_drug_input() -> Optional[str]:
    """Get drug name input from user."""
    try:
        print("\nEnter drug name: ", end='')
        sys.stdout.flush()
        drug_name = input().strip()
        if not drug_name:
            return None
        return drug_name
    except (EOFError, KeyboardInterrupt):
        print("\nInput cancelled")
        return None

def get_dose_input() -> Optional[float]:
    """Get dose input from user."""
    while True:
        try:
            print("Enter dose (mg): ", end='')
            sys.stdout.flush()
            dose_str = input().strip()
            if not dose_str:
                return None
            return float(dose_str)
        except ValueError:
            print("Error: Please enter a valid number for dose")
        except (EOFError, KeyboardInterrupt):
            print("\nInput cancelled")
            return None

def display_drug_info(info: Dict[str, Any]):
    """Display drug information."""
    print("\n=== Drug Properties ===")
    print(f"Molecular Weight: {info.get('molecular_weight', 'N/A')} g/mol")
    print(f"LogP: {info.get('logp', 'N/A')}")
    
    print("\n=== hERG Data ===")
    print(f"IC50: {info.get('herg_ic50', 'N/A')} μM")
    safety_margin = info.get('safety_margin')
    if safety_margin:
        risk = "High Risk" if safety_margin < 30 else "Medium Risk" if safety_margin < 100 else "Low Risk"
        print(f"Safety Margin: {safety_margin:.1f}x ({risk})")

def display_concentration_results(results: Dict[str, float]):
    """Display drug concentration results."""
    print("\n=== Concentration Analysis ===")
    print(f"Theoretical Concentration: {results['theoretical_conc']:.3f} μM")
    print(f"Bioavailable Concentration: {results['bioavailable_conc']:.3f} μM")

def display_literature_results(papers_found: int, qt_points: int, hr_points: int):
    """Display literature analysis results."""
    print("\n=== Literature Analysis ===")
    print(f"Papers Analyzed: {papers_found}")
    print(f"QT Measurements: {qt_points}")
    print(f"HR Measurements: {hr_points}")

def display_loading(text: str):
    """Display loading message."""
    print(f"\n{text}...")
    sys.stdout.flush()

def display_error(text: str):
    """Display error message."""
    print(f"\nError: {text}")

def display_success(text: str):
    """Display success message."""
    print(f"\nSuccess: {text}")

def display_warning(text: str):
    """Display warning message."""
    print(f"\nWarning: {text}")

def display_info(text: str):
    """Display info message."""
    print(f"\nInfo: {text}")
