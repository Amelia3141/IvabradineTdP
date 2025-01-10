"""
Streamlit app for analyzing TdP risk of Ivabradine
"""

import streamlit as st
import plotly.graph_objects as go
import requests
import re
import pandas as pd
import time
import json
import logging
from typing import Dict, List, Optional, Tuple
from ivablib.herg_analyzer import DrugAnalyzer

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config
st.set_page_config(
    page_title="TdP Risk Analyzer",
    page_icon="⚡",
    layout="wide"
)

# Title and description
st.title("⚡ TdP Risk Analyzer")
st.markdown("""
This app analyzes the Torsades de Pointes (TdP) risk of drugs by:
1. Retrieving hERG IC50 values from ChEMBL
2. Checking CredibleMeds database for known risks
3. Searching literature for relevant papers
""")

# Initialize session state
if 'analyzer' not in st.session_state:
    st.session_state.analyzer = DrugAnalyzer(
        email=st.secrets["NCBI_EMAIL"],
        api_key=st.secrets["NCBI_API_KEY"]
    )
if 'drug_name' not in st.session_state:
    st.session_state.drug_name = ""

# Sidebar inputs
with st.sidebar:
    st.header("Analysis Parameters")
    
    # Drug selection
    col1, col2 = st.columns([2,1])
    with col1:
        drug_name = st.text_input("Enter drug name:", value="", key="drug_input")
        analyze_button = st.button("Analyze")

# Main content
if analyze_button and drug_name:
    try:
        with st.spinner("Analyzing drug..."):
            logger.info(f"Starting analysis for {drug_name}")
            analysis = st.session_state.analyzer.analyze_drug(drug_name)
            logger.info(f"Analysis complete: {analysis}")

            if "error" in analysis:
                st.error(f"Error: {analysis['error']}")
            else:
                # Display results in tabs
                tab1, tab2 = st.tabs(["Drug Analysis", "Literature Review"])
                
                with tab1:
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Drug Information")
                        if analysis.get('herg_ic50'):
                            st.write(f"**hERG IC50:** {analysis['herg_ic50']:.2f} μM")
                        else:
                            st.write("**hERG IC50:** Not found")
                        st.write(f"**Source:** {analysis.get('source', 'Unknown')}")
                    
                    with col2:
                        st.subheader("Risk Assessment")
                        if analysis.get('crediblemeds_risk'):
                            risk_text = analysis.get('risk_category', 'Known Risk')
                            st.error(f"⚠️ {risk_text} (CredibleMeds)")
                            st.markdown("[View on CredibleMeds](https://crediblemeds.org)")
                        elif analysis.get('herg_ic50') and analysis['herg_ic50'] < 10:  # Less than 10 μM is concerning
                            st.warning("⚠️ Potential hERG binding detected")
                        else:
                            st.success("✅ No significant hERG binding predicted")
                
                with tab2:
                    st.subheader("Literature Review")
                    literature = analysis.get('literature', pd.DataFrame())
                    if literature.empty:
                        st.info("No relevant literature found.")
                    else:
                        for _, paper in literature.iterrows():
                            with st.expander(paper['title']):
                                st.write(f"**Authors:** {paper['authors']}")
                                st.write(f"**Journal:** {paper['journal']} ({paper['year']})")
                                if paper['abstract']:
                                    st.write(f"**Abstract:** {paper['abstract']}")
                                st.write(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/)")
    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
else:
    st.info("Click 'Analyze' to start the analysis")
