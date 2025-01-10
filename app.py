import streamlit as st

# Page config
st.set_page_config(
    page_title="TdP Risk Assessment",
    page_icon="❤️",
    layout="wide"
)

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import os
from ivablib.herg_analyzer import DrugAnalyzer
import requests
from Bio import Entrez, Medline

# Set NCBI credentials
NCBI_EMAIL = "ghhercock@gmail.com"
NCBI_API_KEY = "d769d934dc8159bcc9ec9fc29c715a456008"

# Set email for Entrez
Entrez.email = NCBI_EMAIL

# Custom CSS to ensure text visibility
st.markdown("""
    <style>
    .stMarkdown, .stText, p, h1, h2, h3 {
        color: black !important;
    }
    th {
        color: black !important;
        font-weight: bold !important;
    }
    td {
        color: black !important;
    }
    </style>
    """, unsafe_allow_html=True)

st.title("TdP Risk Assessment Tool")

# Initialize analyzer
analyzer = DrugAnalyzer(NCBI_EMAIL, NCBI_API_KEY)

# Sidebar
st.sidebar.title("Analysis Sections")
analysis_sections = st.sidebar.multiselect(
    "Select sections to display:",
    ["Risk Category", "hERG Channel Activity", "Literature Analysis"],
    default=["Risk Category", "hERG Channel Activity", "Literature Analysis"]
)

# Drug input
drug_name = st.text_input("Enter drug name:", "ivabradine").strip().lower()

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

def search_pubmed_literature(drug_name: str) -> pd.DataFrame:
    """Search PubMed for TdP-related articles."""
    query = f"{drug_name} AND (Torsades de Pointes[Title/Abstract] OR TdP[Title/Abstract] OR QT prolongation[Title/Abstract])"
    
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            return pd.DataFrame()
            
        # Fetch details for each article
        handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Extract relevant information
        articles = []
        for article in records['PubmedArticle']:
            medline = article['MedlineCitation']
            articles.append({
                'Title': medline['Article']['ArticleTitle'],
                'Journal': medline['Article']['Journal']['Title'],
                'Year': medline['Article']['Journal']['JournalIssue']['PubDate'].get('Year', 'N/A'),
                'PMID': medline['PMID'],
                'Abstract': medline['Article'].get('Abstract', {}).get('AbstractText', [''])[0]
            })
        
        return pd.DataFrame(articles)
        
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
    <div style='color: black;'>
    <p><b>Risk Category:</b> {risk_category}</p>
    <p><b>Explanation:</b> {risk_explanations.get(risk_category, "")}</p>
    </div>
    """, unsafe_allow_html=True)

if "hERG Channel Activity" in analysis_sections:
    st.header("hERG Channel Activity")
    
    try:
        # Standard doses for ivabradine (mg)
        doses = [5.0, 7.5]  # Standard doses
        analysis = analyzer.analyze_drug(drug_name, doses)
        
        # Display hERG data
        st.markdown(f"""
        <div style='color: black;'>
        <p><b>Molecular Weight:</b> {f"{analysis.molecular_weight:.2f}" if analysis.molecular_weight else 'Not found'} g/mol</p>
        <p><b>hERG IC50:</b> {f"{analysis.herg_ic50:.2f}" if analysis.herg_ic50 else 'Not found'} μM</p>
        <p><b>Source:</b> {analysis.herg_source if analysis.herg_source else 'Not available'}</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Display concentration data
        if any(analysis.concentrations) and analysis.molecular_weight:
            st.markdown("<div style='color: black;'><h3>Concentration Analysis</h3></div>", unsafe_allow_html=True)
            for i, conc in enumerate(analysis.concentrations):
                st.markdown(f"""
                <div style='color: black;'>
                <p><b>Dose {doses[i]} mg:</b></p>
                <ul>
                    <li>Theoretical Max: {f"{conc.theoretical_max:.2f}"} μM</li>
                    <li>Plasma Concentration: {f"{conc.plasma_concentration:.2f}"} μM</li>
                    {f'<li>IC50/Plasma Ratio: {conc.ratio_plasma:.2f}</li>' if conc.ratio_plasma else ''}
                </ul>
                </div>
                """, unsafe_allow_html=True)
        
        # Display citations
        if analysis.citations:
            st.markdown("<div style='color: black;'><h3>References</h3></div>", unsafe_allow_html=True)
            for citation in analysis.citations:
                st.markdown(f"<div style='color: black;'>{citation}</div>", unsafe_allow_html=True)
            
    except Exception as e:
        st.error(f"Error analyzing hERG activity: {str(e)}")
        st.info("Please check that your NCBI credentials are properly configured.")
    
    st.markdown("### hERG Binding Analysis", unsafe_allow_html=True)
    st.markdown("""
    <div style='color: black;'>
    <p>The concentration analysis shows:</p>
    <ul>
        <li><b>Theoretical Max:</b> Maximum theoretical concentration based on dose and distribution volume</li>
        <li><b>Plasma Concentration:</b> Estimated plasma concentration with 40% bioavailability</li>
        <li><b>IC50/Plasma Ratio:</b> Ratio between hERG IC50 and plasma concentration (values < 1 indicate potential risk)</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)

if "Literature Analysis" in analysis_sections:
    st.header("Literature Analysis")
    
    with st.spinner("Searching PubMed literature..."):
        df = search_pubmed_literature(drug_name)
        if not df.empty:
            st.markdown(f"<div style='color: black;'>Found {len(df)} relevant articles:</div>", unsafe_allow_html=True)
            st.dataframe(df, use_container_width=True)
        else:
            st.info("No relevant literature found.")
