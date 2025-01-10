from ivabapp2 import plot_combined_analysis

# Sample diphenhydramine case reports
papers = [
    {
        'title': 'Diphenhydramine-induced QT prolongation and torsade de pointes: An uncommon effect of a common drug',
        'abstract': '''A 67-year-old female presented with palpitations after taking diphenhydramine 50mg. 
        Initial QTc was 510ms with HR 72bpm. Blood pressure was 125/80 mmHg. 
        Patient had history of hypertension and was on amlodipine. 
        Treatment involved discontinuation of diphenhydramine and magnesium supplementation.''',
        'drug_concentration': '0.15',  # μM
        'herg_ic50': '0.5',  # μM
        'naranjo_score': 7,
        'tisdale_score': 8,
        'who_umc': 'Probable'
    },
    {
        'title': 'Torsades de pointes associated with high-dose diphenhydramine ingestion',
        'abstract': '''A 42-year-old male presented after intentional ingestion of diphenhydramine (estimated 800mg). 
        QTc prolonged to 580ms from baseline 440ms. Heart rate dropped to 45bpm. BP 110/70 mmHg. 
        Patient developed TdP requiring cardioversion. Treatment included activated charcoal and magnesium.''',
        'drug_concentration': '2.4',  # μM
        'herg_ic50': '0.5',  # μM
        'naranjo_score': 9,
        'tisdale_score': 11,
        'who_umc': 'Certain'
    }
]

# QT and HR data points for nomogram
qt_data = [
    {'qt': 510, 'hr': 72},
    {'qt': 580, 'hr': 45}
]

# Generate analysis
plot_combined_analysis(papers, "Diphenhydramine", qt_data)
