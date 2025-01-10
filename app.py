"""Streamlit app for analyzing TdP risk."""
import os
import sys
from typing import Dict, List, Any, Optional
import streamlit as st
import plotly.graph_objects as go
import requests
from Bio import Entrez, Medline

# Add the project root directory to Python path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.append(project_root)

from ivablib.herg_analyzer import DrugAnalyzer
from ivablib.case_report_analyzer import CaseReportAnalyzer

# Set NCBI credentials
Entrez.email = "amelia.glasauer@gmail.com"
Entrez.api_key = "8b23b0c60e6b8c8a40f2867c88a9d43c5f09"

# Page config
st.set_page_config(
    page_title="TdP Risk Assessment",
    page_icon="❤️",
    layout="wide"
)

# Custom CSS to ensure text visibility
st.markdown("""
    <style>
    .stMarkdown, .stText, p, h1, h2, h3, h4 {
        color: white !important;
    }
    .stDataFrame {
        color: white !important;
    }
    .stDataFrame td, .stDataFrame th {
        color: white !important;
        background-color: rgba(255, 255, 255, 0.1) !important;
    }
    .stDataFrame th {
        font-weight: bold !important;
    }
    .risk-info {
        background-color: rgba(255, 255, 255, 0.1);
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
    }
    .herg-data {
        background-color: rgba(255, 255, 255, 0.1);
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
    }
    .concentration-data {
        background-color: rgba(255, 255, 255, 0.1);
        padding: 15px;
        border-radius: 8px;
        margin: 5px 0;
    }
    .literature-data {
        background-color: rgba(255, 255, 255, 0.1);
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
    }
    </style>
    """, unsafe_allow_html=True)

# Initialize analyzers
drug_analyzer = DrugAnalyzer(Entrez.email, Entrez.api_key)

# Main content
st.title("TdP Risk Assessment")
st.markdown("Analyze the risk of Torsades de Pointes (TdP) for a given drug.")

# Drug input
drug_name = st.text_input("Enter drug name:", "Ivabradine").strip().lower()

# Analysis sections selection
analysis_sections = st.multiselect(
    "Select analysis sections:",
    ["Risk Category", "hERG Channel Activity", "Literature Analysis"],
    default=["Risk Category", "hERG Channel Activity", "Literature Analysis"]
)

def check_drug_risk(drug_name):
    """Check drug risk category in CredibleMeds database."""
    risk_categories = {
        "ivabradine": "Known Risk of TdP",
        "amiodarone": "Known Risk of TdP",
        "sotalol": "Known Risk of TdP",
        "methadone": "Known Risk of TdP",
        "donepezil": "Conditional Risk of TdP",
        "galantamine": "Conditional Risk of TdP",
        "rivastigmine": "No Known Risk of TdP"
    }
    return risk_categories.get(drug_name, "Not Found in Database")

def extract_clinical_data(text: str) -> Dict:
    """Extract clinical data from text using regex patterns."""
    data = {}
    
    # Age extraction
    age_pattern = r'(\d+)[\s-]*(year|yr|y)[s\s-]*old'
    age_match = re.search(age_pattern, text.lower())
    data['Age'] = age_match.group(1) if age_match else None
    
    # Sex extraction
    sex_pattern = r'\b(male|female)\b'
    sex_match = re.search(sex_pattern, text.lower())
    data['Sex'] = sex_match.group(1).capitalize() if sex_match else None
    
    # Dose extraction
    dose_pattern = r'(\d+(?:\.\d+)?)\s*mg'
    dose_match = re.search(dose_pattern, text.lower())
    data['Oral Dose (mg)'] = dose_match.group(1) if dose_match else None
    
    # QT interval extraction
    qt_pattern = r'QT(?:c)?\s*(?:interval)?\s*(?:of|was|is|=)?\s*(\d+)(?:\s*±\s*\d+)?\s*(?:ms|msec)'
    qt_match = re.search(qt_pattern, text.lower())
    data['Uncorrected QT (ms)'] = qt_match.group(1) if qt_match else None
    
    # QTc extraction
    qtc_pattern = r'QTc\s*(?:interval)?\s*(?:of|was|is|=)?\s*(\d+)(?:\s*±\s*\d+)?\s*(?:ms|msec)'
    qtc_match = re.search(qtc_pattern, text.lower())
    data['QTc'] = qtc_match.group(1) if qtc_match else None
    
    # Heart rate extraction
    hr_pattern = r'(?:heart rate|hr|pulse)\s*(?:of|was|is|=)?\s*(\d+)(?:\s*±\s*\d+)?\s*(?:bpm|beats per minute)'
    hr_match = re.search(hr_pattern, text.lower())
    data['Heart Rate (bpm)'] = hr_match.group(1) if hr_match else None
    
    # Blood pressure extraction
    bp_pattern = r'(?:blood pressure|bp)\s*(?:of|was|is|=)?\s*(\d+/\d+)(?:\s*±\s*\d+)?\s*(?:mmHg)'
    bp_match = re.search(bp_pattern, text.lower())
    data['Blood Pressure (mmHg)'] = bp_match.group(1) if bp_match else None
    
    # TdP presence
    data['Torsades de Pointes?'] = 'Yes' if re.search(r'\b(?:torsade|tdp)\b', text.lower()) else 'No'
    
    return data

def search_pubmed_literature(drug_name: str) -> pd.DataFrame:
    """Search PubMed for TdP-related articles."""
    query = f"{drug_name} AND (Torsades de Pointes[Title/Abstract] OR TdP[Title/Abstract] OR QT prolongation[Title/Abstract])"
    
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)  
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            return pd.DataFrame()
            
        # Fetch details for each article
        handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Extract relevant information
        papers = []
        for article in records['PubmedArticle']:
            try:
                medline = article['MedlineCitation']
                article_data = medline['Article']
                
                # Extract authors
                authors = article_data.get('AuthorList', [])
                author_names = []
                for author in authors:
                    if isinstance(author, dict):
                        last_name = author.get('LastName', '')
                        fore_name = author.get('ForeName', '')
                        author_names.append(f"{last_name} {fore_name}")
                
                # Get publication date
                pub_date = article_data['Journal']['JournalIssue']['PubDate']
                year = pub_date.get('Year', 'N/A')
                
                # Extract abstract text safely
                abstract = ''
                abstract_data = article_data.get('Abstract', {})
                if abstract_data:
                    abstract_text = abstract_data.get('AbstractText', [])
                    if abstract_text:
                        if isinstance(abstract_text, list):
                            # Handle labeled sections
                            sections = []
                            for section in abstract_text:
                                if isinstance(section, str):
                                    sections.append(section)
                                elif hasattr(section, 'attributes') and hasattr(section, 'string'):
                                    label = section.attributes.get('Label', '')
                                    text = str(section)
                                    sections.append(f"{label}: {text}" if label else text)
                            abstract = ' '.join(sections)
                        else:
                            abstract = str(abstract_text)
                
                papers.append({
                    'title': article_data['ArticleTitle'],
                    'authors': ', '.join(author_names) if author_names else 'No authors listed',
                    'year': year,
                    'journal': article_data['Journal'].get('Title', 'Journal not specified'),
                    'pmid': medline['PMID'],
                    'abstract': abstract
                })
            except Exception as e:
                st.warning(f"Error processing article: {str(e)}")
                continue
        
        if not papers:
            return pd.DataFrame()
            
        # Initialize case report analyzer
        analyzer = CaseReportAnalyzer()
        
        # Analyze papers and create DataFrame
        df = analyzer.analyze_papers(papers, drug_name)
        
        # Reorder columns to match desired order
        column_order = [
            'title', 'age', 'sex', 'qtc', 'qt_uncorrected',
            'heart_rate', 'blood_pressure', 'had_tdp',
            'patient_type', 'treatment_successful', 'treatment_duration',
            'drug_combinations', 'authors', 'year', 'journal', 'pmid'
        ]
        
        # Only include columns that exist in the DataFrame
        existing_columns = [col for col in column_order if col in df.columns]
        df = df[existing_columns]
        
        # Rename columns for display
        column_names = {
            'title': 'Case Report Title',
            'age': 'Age',
            'sex': 'Sex',
            'qtc': 'QTc (ms)',
            'qt_uncorrected': 'Uncorrected QT (ms)',
            'heart_rate': 'Heart Rate (bpm)',
            'blood_pressure': 'Blood Pressure (mmHg)',
            'had_tdp': 'Torsades de Pointes?',
            'patient_type': 'Patient Type',
            'treatment_successful': 'Treatment Successful',
            'treatment_duration': 'Treatment Duration',
            'drug_combinations': 'Drug Combinations',
            'authors': 'Authors',
            'year': 'Year',
            'journal': 'Journal',
            'pmid': 'PMID'
        }
        df = df.rename(columns=column_names)
        
        return df
        
    except Exception as e:
        st.error(f"Error searching PubMed: {str(e)}")
        return pd.DataFrame()

def create_risk_gauge(risk_category):
    """Create a Plotly gauge chart for risk visualization."""
    risk_levels = {
        "Known Risk of TdP": 1.0,
        "Possible Risk of TdP": 0.66,
        "Conditional Risk of TdP": 0.33,
        "No Known Risk of TdP": 0.0,
        "Not Found in Database": 0.0
    }
    
    value = risk_levels.get(risk_category, 0)
    
    fig = go.Figure(go.Indicator(
        mode = "gauge",
        value = value,
        gauge = {
            'axis': {'visible': False, 'range': [0, 1]},
            'bar': {'color': "rgba(0,0,0,0)"},
            'steps': [
                {'range': [0, 0.33], 'color': "lightgreen"},
                {'range': [0.33, 0.66], 'color': "yellow"},
                {'range': [0.66, 1], 'color': "red"}
            ],
            'threshold': {
                'line': {'color': "black", 'width': 4},
                'thickness': 0.75,
                'value': value
            }
        }
    ))
    
    fig.update_layout(
        title={
            'text': f"TdP Risk<br><a href='https://crediblemeds.org/' target='_blank' style='color:blue;'>Data from CredibleMeds</a>",
            'y':0.8,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        height=300,
        margin=dict(l=30, r=30, t=100, b=30)
    )
    return fig

def analyze_herg_activity(drug_name: str, doses: List[float]) -> Dict[str, Any]:
    """Analyze hERG channel activity for drug."""
    try:
        # Initialize analyzer with NCBI credentials
        analyzer = DrugAnalyzer(
            email=Entrez.email,
            api_key=Entrez.api_key
        )
        
        # Analyze drug
        analysis = analyzer.analyze_drug(drug_name, doses)
        if not analysis:
            return {
                "error": "No hERG IC50 data found in literature.",
                "concentrations": []
            }
            
        # Format concentrations for display
        concentrations = []
        for i, conc in enumerate(analysis.concentrations):
            conc_data = {
                "dose": doses[i],
                "theoretical_max": f"{conc.theoretical_max:.2f}",
                "plasma_concentration": f"{conc.plasma_concentration:.2f}",
                "ratio_theoretical": f"{conc.ratio_theoretical:.2f}" if conc.ratio_theoretical else "N/A",
                "ratio_plasma": f"{conc.ratio_plasma:.2f}" if conc.ratio_plasma else "N/A"
            }
            concentrations.append(conc_data)
            
        return {
            "name": analysis.name,
            "molecular_weight": f"{analysis.molecular_weight:.2f}" if analysis.molecular_weight else "N/A",
            "herg_ic50": f"{analysis.herg_ic50:.2f}" if analysis.herg_ic50 else "N/A",
            "herg_source": analysis.herg_source or "N/A",
            "theoretical_binding": analysis.theoretical_binding,
            "concentrations": concentrations,
            "citations": analysis.citations
        }
        
    except Exception as e:
        return {
            "error": str(e),
            "concentrations": []
        }

# Main content
if "Risk Category" in analysis_sections:
    st.header("Risk Category Analysis")
    risk_category = check_drug_risk(drug_name)
    st.plotly_chart(create_risk_gauge(risk_category), use_container_width=True)
    
    # Display risk category explanation
    risk_explanations = {
        "Known Risk of TdP": "Substantial evidence supports the conclusion that these drugs prolong QT intervals and have a risk of TdP when used as directed in labeling.",
        "Possible Risk of TdP": "Substantial evidence supports the conclusion that these drugs can cause QT prolongation but there is insufficient evidence that they, when used as directed in labeling, have a risk of causing TdP.",
        "Conditional Risk of TdP": "Substantial evidence supports the conclusion that these drugs prolong QT and have a risk of developing TdP but only under certain known conditions.",
        "No Known Risk of TdP": "Available evidence does not support the conclusion that these drugs prolong QT intervals nor have a risk of TdP when used as directed in labeling.",
        "Not Found in Database": "This drug was not found in the CredibleMeds database. Please verify the drug name or consult additional sources."
    }
    
    st.markdown(f"""
    <div class='risk-info'>
    <p><b>Risk Category:</b> {risk_category}</p>
    <p><b>Explanation:</b> {risk_explanations.get(risk_category, "")}</p>
    </div>
    """, unsafe_allow_html=True)

if "hERG Channel Activity" in analysis_sections:
    st.header("hERG Channel Activity")
    
    try:
        # Standard doses for ivabradine (mg)
        doses = [5.0, 7.5]  # Standard doses
        analysis = analyze_herg_activity(drug_name, doses)
        
        if analysis.get("error"):
            st.error(analysis["error"])
        else:
            # Display hERG data
            st.markdown("<div class='herg-data'>", unsafe_allow_html=True)
            
            if analysis["herg_ic50"] != "N/A":
                st.markdown(f"**hERG IC50:** {analysis['herg_ic50']}")
                if analysis["herg_source"] != "N/A":
                    st.markdown(f"**Source:** {analysis['herg_source']}")
            else:
                st.warning("No hERG IC50 data found in literature.")
            
            st.markdown("</div>", unsafe_allow_html=True)
            
            # Display concentration analysis
            st.subheader("hERG Binding Analysis")
            st.markdown("<div class='concentration-data'>", unsafe_allow_html=True)
            st.markdown("The concentration analysis shows:")
            
            for conc in analysis["concentrations"]:
                st.markdown(f"""
                * **Dose:** {conc['dose']} mg
                * **Theoretical Max:** {conc['theoretical_max']} μM (Maximum theoretical concentration based on dose and distribution volume)
                * **Plasma Concentration:** {conc['plasma_concentration']} μM (Estimated plasma concentration with 40% bioavailability)
                * **IC50/Theoretical Ratio:** {conc['ratio_theoretical']} (Ratio between hERG IC50 and theoretical concentration - values < 1 indicate potential risk)
                * **IC50/Plasma Ratio:** {conc['ratio_plasma']} (Ratio between hERG IC50 and plasma concentration - values < 1 indicate potential risk)
                """)
            
            st.markdown("</div>", unsafe_allow_html=True)
            
            # Citations
            if analysis["citations"]:
                st.markdown("<div class='citations'>", unsafe_allow_html=True)
                st.markdown("**References:**")
                for citation in analysis["citations"]:
                    st.markdown(f"* {citation}")
                st.markdown("</div>", unsafe_allow_html=True)
            
    except Exception as e:
        st.error(f"Error analyzing hERG activity: {str(e)}")
        st.info("Please check that your NCBI credentials are properly configured.")

if "Literature Analysis" in analysis_sections:
    st.header("Literature Analysis")
    
    with st.spinner("Searching PubMed literature..."):
        df = search_pubmed_literature(drug_name)
        if not df.empty:
            st.markdown(f"<div style='color: white;'>Found {len(df)} relevant case reports:</div>", unsafe_allow_html=True)
            st.markdown("<div class='literature-data'>", unsafe_allow_html=True)
            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
            )
            st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.info("No relevant literature found.")
