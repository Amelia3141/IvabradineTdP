from ivabapp2 import plot_combined_analysis

# Sample paper data for ivabradine
papers = [
    {
        'title': 'Ivabradine-induced Torsades de Pointes in a Patient with Previous Myocardial Infarction',
        'year': '2023',
        'authors': 'Cardio Research Team',
        'abstract': 'A 72-year-old male with history of MI presented with QTc prolongation of 510ms and heart rate of 40 bpm while on ivabradine therapy. Initial QTc was 440ms with HR 75 bpm. Blood pressure was 125/80 mmHg.',
        'naranjo_score': 8,
        'tisdale_score': 11,
        'who_umc': 'Probable'
    },
    {
        'title': 'QT Interval Changes with Ivabradine Treatment',
        'year': '2022',
        'authors': 'Heart Rhythm Study Group',
        'abstract': 'Patient developed QTc of 485ms with heart rate 50 bpm during ivabradine treatment. Baseline QTc was 430ms with HR 85 bpm.',
        'naranjo_score': 7,
        'tisdale_score': 9,
        'who_umc': 'Possible'
    },
    {
        'title': 'Bradycardia and QT Prolongation with Ivabradine',
        'year': '2021',
        'authors': 'Cardiac Safety Team',
        'abstract': 'Case series of three patients showing QTc prolongation (495ms, 475ms, 460ms) with corresponding heart rates of 45, 55, and 60 bpm respectively during ivabradine therapy.',
        'naranjo_score': 6,
        'tisdale_score': 8,
        'who_umc': 'Possible'
    }
]

# QT and HR data points
qt_data = [
    {'qt': 510, 'hr': 40},  # First case initial
    {'qt': 440, 'hr': 75},  # First case baseline
    {'qt': 485, 'hr': 50},  # Second case
    {'qt': 430, 'hr': 85},  # Second case baseline
    {'qt': 495, 'hr': 45},  # Third case patient 1
    {'qt': 475, 'hr': 55},  # Third case patient 2
    {'qt': 460, 'hr': 60}   # Third case patient 3
]

# Generate plot
plot_combined_analysis(papers, "Ivabradine", qt_data)
