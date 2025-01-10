from ivabapp2 import plot_combined_analysis

# Comprehensive paper data for amiodarone
papers = [
    {
        'title': 'Amiodarone-induced Torsade de Pointes in Elderly Patient',
        'year': '2023',
        'authors': 'Cardiac Arrhythmia Research Group',
        'abstract': '''75-year-old female developed TdP after receiving IV amiodarone for atrial fibrillation. 
        Initial QTc was 440ms with HR 98bpm. After loading dose, QTc prolonged to 580ms with HR dropping to 42bpm. 
        Blood pressure was 110/70 mmHg. Patient had hypokalemia (K+ 3.2) and was on concurrent ciprofloxacin.''',
        'naranjo_score': 9,
        'tisdale_score': 14,
        'who_umc': 'Definite',
        'max_qt': 580,
        'tdp_observed': True,
        'drug_concentration': '2.1',
        'herg_ic50': '0.3'
    },
    {
        'title': 'QT Prolongation in Chronic Amiodarone Therapy',
        'year': '2022',
        'authors': 'Electrophysiology Study Team',
        'abstract': '''Long-term follow-up of 3 patients on oral amiodarone maintenance therapy.
        Case 1: QTc 520ms, HR 55bpm
        Case 2: QTc 495ms, HR 48bpm
        Case 3: QTc 510ms, HR 50bpm
        All patients maintained stable rhythm without TdP despite QT prolongation.''',
        'naranjo_score': 7,
        'tisdale_score': 11,
        'who_umc': 'Probable',
        'max_qt': 520,
        'tdp_observed': False,
        'drug_concentration': '1.8',
        'herg_ic50': '0.3'
    },
    {
        'title': 'Acute Amiodarone Loading and QT Dynamics',
        'year': '2021',
        'authors': 'Cardiology Research Unit',
        'abstract': '''Prospective study of 5 patients receiving IV amiodarone loading.
        Patient 1: Baseline QTc 430ms/HR 90bpm → Peak QTc 515ms/HR 45bpm
        Patient 2: Baseline QTc 445ms/HR 85bpm → Peak QTc 490ms/HR 52bpm
        Patient 3: Baseline QTc 420ms/HR 95bpm → Peak QTc 505ms/HR 48bpm
        Patient 4: Baseline QTc 450ms/HR 88bpm → Peak QTc 485ms/HR 58bpm
        Patient 5: Baseline QTc 425ms/HR 92bpm → Peak QTc 495ms/HR 50bpm''',
        'naranjo_score': 8,
        'tisdale_score': 10,
        'who_umc': 'Probable',
        'max_qt': 515,
        'tdp_observed': False,
        'drug_concentration': '1.9',
        'herg_ic50': '0.3'
    }
]

# Comprehensive QT/HR data points including baseline and peak values
qt_data = [
    # First paper
    {'qt': 440, 'hr': 98},  # Initial
    {'qt': 580, 'hr': 42},  # After loading
    
    # Second paper
    {'qt': 520, 'hr': 55},  # Case 1
    {'qt': 495, 'hr': 48},  # Case 2
    {'qt': 510, 'hr': 50},  # Case 3
    
    # Third paper
    {'qt': 430, 'hr': 90},  # Patient 1 baseline
    {'qt': 515, 'hr': 45},  # Patient 1 peak
    {'qt': 445, 'hr': 85},  # Patient 2 baseline
    {'qt': 490, 'hr': 52},  # Patient 2 peak
    {'qt': 420, 'hr': 95},  # Patient 3 baseline
    {'qt': 505, 'hr': 48},  # Patient 3 peak
    {'qt': 450, 'hr': 88},  # Patient 4 baseline
    {'qt': 485, 'hr': 58},  # Patient 4 peak
    {'qt': 425, 'hr': 92},  # Patient 5 baseline
    {'qt': 495, 'hr': 50}   # Patient 5 peak
]

# Generate plot
plot_combined_analysis(papers, "Amiodarone", qt_data)
