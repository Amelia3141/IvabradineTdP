import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import os
import logging
import requests
import urllib3
from bs4 import BeautifulSoup
from Bio import Entrez, Medline
from typing import List, Dict, Optional
import sys
import time

# Set email for Entrez
Entrez.email = "your.email@example.com"

# Page config
st.set_page_config(
    page_title="TdP Risk Assessment",
    page_icon="❤️",
    layout="wide"
)

# Custom CSS
st.markdown("""
    <style>
    .main {
        background-color: #f5f5f5;
    }
    .stButton>button {
        width: 100%;
        background-color: #ff4b4b;
        color: white;
    }
    .stProgress .st-bo {
        background-color: #ff4b4b;
    }
    </style>
    """, unsafe_allow_html=True)

# Title and description
st.title("TdP Risk Assessment Dashboard")
st.markdown("""
This dashboard provides a comprehensive analysis of Torsades de Pointes (TdP) risk for drugs, 
with a focus on Ivabradine. It combines data from multiple sources including CredibleMeds, 
PubMed literature analysis, and hERG channel activity.
""")

# Sidebar
with st.sidebar:
    st.header("Analysis Options")
    drug_name = st.text_input("Drug Name:", value="Ivabradine")
    analysis_sections = st.multiselect(
        "Select Analysis Sections:",
        ["Risk Category", "Literature Analysis", "hERG Activity"],
        default=["Risk Category", "Literature Analysis", "hERG Activity"]
    )
    st.divider()
    st.markdown("### About")
    st.markdown("""
    This tool assesses TdP risk by analyzing:
    - CredibleMeds risk categorization
    - PubMed literature review
    - hERG channel interactions
    """)

# Risk Categories
RISK_CATEGORIES = {
    "Known Risk of TdP": {
        "color": "#FF4444",
        "description": "These drugs prolong the QT interval AND are clearly associated with a known risk of TdP, even when taken as recommended."
    },
    "Possible Risk of TdP": {
        "color": "#FFA500",
        "description": "These drugs can cause QT prolongation BUT there is a lack of evidence for a risk of TdP when taken as recommended."
    },
    "Conditional Risk of TdP": {
        "color": "#FFD700",
        "description": "These drugs are associated with TdP BUT only under certain conditions."
    },
    "Not Listed": {
        "color": "#CCCCCC",
        "description": "This drug is not listed in the CredibleMeds database."
    }
}

def check_drug_risk(drug_name):
    """Check drug risk category in CredibleMeds database."""
    with open("crediblemeds_data.txt", "r") as f:
        data_lines = [line for line in f.readlines() if not line.startswith("#") and line.strip()]
    
    drug_name = drug_name.lower()
    for line in data_lines:
        name, category = [x.strip() for x in line.split("|")]
        if name.lower() == drug_name:
            return category
    return "Not Listed"

def create_risk_gauge(risk_category):
    """Create a Plotly gauge chart for risk visualization."""
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        gauge={
            'axis': {'range': [None, 3], 
                    'ticktext': ['Low', 'Medium', 'High'],
                    'tickvals': [0, 1.5, 3]},
            'bar': {'color': RISK_CATEGORIES[risk_category]["color"]},
            'bgcolor': "white",
            'borderwidth': 2,
            'bordercolor': "gray",
            'steps': [
                {'range': [0, 1], 'color': "#FFD700"},
                {'range': [1, 2], 'color': "#FFA500"},
                {'range': [2, 3], 'color': "#FF4444"}
            ],
            'threshold': {
                'line': {'color': "black", 'width': 4},
                'thickness': 0.75,
                'value': 2.5 if risk_category == "Known Risk of TdP" else 
                         1.5 if risk_category == "Possible Risk of TdP" else
                         0.5 if risk_category == "Conditional Risk of TdP" else 0
            }
        },
        value= 3 if risk_category == "Known Risk of TdP" else 
               2 if risk_category == "Possible Risk of TdP" else
               1 if risk_category == "Conditional Risk of TdP" else 0,
        title={'text': risk_category, 'font': {'size': 24}}
    ))
    
    fig.update_layout(
        height=300,
        margin=dict(t=25, b=25, l=25, r=25)
    )
    return fig

def search_pubmed_case_reports(drug_name):
    """Search PubMed for case reports about the drug."""
    query = f"{drug_name} AND (case report[Title/Abstract] OR case study[Title/Abstract]) AND (QT prolongation[Title/Abstract] OR Torsades de pointes[Title/Abstract] OR TdP[Title/Abstract])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record['IdList']:
        return pd.DataFrame()
        
    handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
    records = Medline.parse(handle)
    return pd.DataFrame(list(records))

def search_pubmed_cohort_studies(drug_name):
    """Search PubMed for cohort studies about the drug."""
    query = f"{drug_name} AND (cohort study[Title/Abstract]) AND (QT prolongation[Title/Abstract] OR Torsades de pointes[Title/Abstract] OR TdP[Title/Abstract])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record['IdList']:
        return pd.DataFrame()
        
    handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
    records = Medline.parse(handle)
    return pd.DataFrame(list(records))

def search_pubmed_clinical_trials(drug_name):
    """Search PubMed for clinical trials about the drug."""
    query = f"{drug_name} AND (clinical trial[Title/Abstract]) AND (QT prolongation[Title/Abstract] OR Torsades de pointes[Title/Abstract] OR TdP[Title/Abstract])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    if not record['IdList']:
        return pd.DataFrame()
        
    handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
    records = Medline.parse(handle)
    return pd.DataFrame(list(records))

def get_herg_activity(drug_name):
    """Get hERG activity data for the drug."""
    # Add implementation for hERG activity data retrieval
    pass

# Main content
if "Risk Category" in analysis_sections:
    st.header("Risk Category Analysis")
    risk_category = check_drug_risk(drug_name)
    
    col1, col2 = st.columns([1, 1])
    with col1:
        st.plotly_chart(create_risk_gauge(risk_category), use_container_width=True)
    with col2:
        st.markdown("### Risk Description")
        st.info(RISK_CATEGORIES[risk_category]["description"])
        st.markdown("### Risk Factors")
        st.markdown("""
        - QT interval prolongation
        - Electrolyte imbalances
        - Underlying heart conditions
        - Drug interactions
        """)

if "Literature Analysis" in analysis_sections:
    st.header("Literature Analysis")
    
    tabs = st.tabs(["Case Reports", "Cohort Studies", "Clinical Trials"])
    
    with tabs[0]:
        with st.spinner("Analyzing case reports..."):
            case_reports = search_pubmed_case_reports(drug_name)
            st.dataframe(case_reports, use_container_width=True)
    
    with tabs[1]:
        with st.spinner("Analyzing cohort studies..."):
            cohort_studies = search_pubmed_cohort_studies(drug_name)
            st.dataframe(cohort_studies, use_container_width=True)
    
    with tabs[2]:
        with st.spinner("Analyzing clinical trials..."):
            clinical_trials = search_pubmed_clinical_trials(drug_name)
            st.dataframe(clinical_trials, use_container_width=True)

if "hERG Activity" in analysis_sections:
    st.header("hERG Channel Activity")
    with st.spinner("Analyzing hERG activity..."):
        try:
            herg_data = get_herg_activity(drug_name)
            st.code(herg_data, language="text")
            
            # Add visualization placeholder for hERG data
            st.markdown("### hERG Binding Analysis")
            st.info("Visualization of hERG channel binding characteristics will be added in future updates.")
        except Exception as e:
            st.error("Error analyzing hERG activity: " + str(e))
