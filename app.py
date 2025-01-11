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
from ivablib.herg_analyzer import DrugAnalyzer
from ivablib.pubmed4125 import analyze_literature

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config
st.set_page_config(
    page_title="TdP Risk Analyzer",
    page_icon="⚡",
    layout="wide"
)

# Initialize session state
if 'analyzer' not in st.session_state:
    st.session_state.analyzer = DrugAnalyzer(
        email=st.secrets["NCBI_EMAIL"],
        api_key=st.secrets["NCBI_API_KEY"]
    )

# Title
st.title("⚡ TdP Risk Analyzer")
st.markdown("""
This app analyzes the Torsades de Pointes (TdP) risk of drugs by:
1. Retrieving hERG IC50 values from ChEMBL
2. Checking CredibleMeds database for known risks
3. Searching literature for relevant papers
""")

# Input
col1, col2 = st.columns([2,1])
with col1:
    drug_name = st.text_input("Enter drug name:", value="", key="drug_input")
    analyze_button = st.button("Analyze")

# Main content
if analyze_button:
    try:
        with st.spinner("Analyzing drug..."):
            logger.info(f"Starting analysis for {drug_name}")
            
            # Get drug analysis
            analysis = st.session_state.analyzer.analyze_drug(drug_name)
            logger.info(f"Drug analysis complete")
            
            # Get literature analysis
            literature = analyze_literature(drug_name)
            analysis['literature'] = literature
            logger.info(f"Literature analysis complete")

            if "error" in analysis:
                st.error(f"Error: {analysis['error']}")
            else:
                # Display results in tabs
                tab1, tab2 = st.tabs(["Drug Analysis", "Literature Review"])
                
                with tab1:
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Drug Information")
                        st.write(f"**PubChem CID:** {analysis.get('cid', 'Not found')}")
                        st.write(f"**Molecular Weight:** {analysis.get('molecular_weight', 'Not found')} g/mol")
                        st.write(f"**hERG IC50:** {analysis.get('herg_ic50', 'None')} μM")
                        st.write(f"**Source:** {analysis.get('source', 'Unknown')}")
                        
                        # Show concentration analysis
                        concentrations = analysis.get('concentrations')
                        if concentrations:
                            st.subheader("Concentration Analysis")
                            st.write(f"**Dose:** {concentrations['dose_mg']} mg")
                            st.write(f"**Theoretical Max:** {concentrations['theoretical_uM']:.2f} μM")
                            st.write(f"**Plasma Concentration:** {concentrations['plasma_uM']:.2f} μM")
                            
                            if concentrations['theoretical_ratio']:
                                st.write(f"**Theoretical/IC50 Ratio:** {concentrations['theoretical_ratio']:.2f}")
                                if concentrations['theoretical_ratio'] > 1:
                                    st.warning("⚠️ Theoretical concentration exceeds hERG IC50")
                                elif concentrations['theoretical_ratio'] > 0.1:
                                    st.warning("⚠️ Theoretical concentration approaches hERG IC50")
                            
                            if concentrations['plasma_ratio']:
                                st.write(f"**Plasma/IC50 Ratio:** {concentrations['plasma_ratio']:.2f}")
                                if concentrations['plasma_ratio'] > 1:
                                    st.error("⚠️ Plasma concentration exceeds hERG IC50")
                                elif concentrations['plasma_ratio'] > 0.1:
                                    st.warning("⚠️ Plasma concentration approaches hERG IC50")
                    
                    with col2:
                        st.subheader("Risk Assessment")
                        if analysis.get('crediblemeds_risk'):
                            risk_text = analysis.get('risk_category', 'Known Risk')
                            st.error(f"⚠️ {risk_text} (CredibleMeds)")
                            st.markdown("[View on CredibleMeds](https://crediblemeds.org)")
                        elif analysis.get('theoretical_binding'):
                            st.warning("⚠️ Potential hERG binding detected")
                            if concentrations and concentrations['theoretical_ratio']:
                                st.write(f"Theoretical concentration is {concentrations['theoretical_ratio']:.1f}x the IC50")
                        else:
                            st.success("✅ No significant hERG binding predicted")
                
                with tab2:
                    st.subheader("Literature Review")
                    literature = analysis.get('literature', {})
                    
                    if 'error' in literature:
                        st.error(f"Error analyzing literature: {literature['error']}")
                    else:
                        # Show summary
                        summary = literature.get('summary', {})
                        st.write(f"Found {summary.get('total_papers', 0)} papers ({summary.get('papers_with_text', 0)} with full text)")
                        
                        # Show key findings
                        if summary.get('key_findings'):
                            st.subheader("Key Findings")
                            for finding in summary['key_findings']:
                                st.write(f"- {finding}")
                        
                        # Show papers by type
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.subheader("Case Reports")
                            for paper in literature.get('case_reports', []):
                                with st.expander(paper['Title']):
                                    st.write(f"**Authors:** {paper['Authors']}")
                                    st.write(f"**Journal:** {paper['Journal']} ({paper['Year']})")
                                    st.write(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper['PMID']}/)")
                        
                        with col2:
                            st.subheader("Cohort Studies")
                            for paper in literature.get('cohort_studies', []):
                                with st.expander(paper['Title']):
                                    st.write(f"**Authors:** {paper['Authors']}")
                                    st.write(f"**Journal:** {paper['Journal']} ({paper['Year']})")
                                    st.write(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper['PMID']}/)")
                        
                        with col3:
                            st.subheader("Clinical Trials")
                            for paper in literature.get('clinical_trials', []):
                                with st.expander(paper['Title']):
                                    st.write(f"**Authors:** {paper['Authors']}")
                                    st.write(f"**Journal:** {paper['Journal']} ({paper['Year']})")
                                    st.write(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper['PMID']}/)")
                        
                        # Show full text findings
                        if literature.get('full_texts'):
                            st.subheader("Full Text Analysis")
                            for paper in literature['full_texts']:
                                with st.expander(paper['title']):
                                    st.write(f"**Authors:** {', '.join(paper['authors'])}")
                                    st.write(f"**Journal:** {paper['journal']} ({paper['year']})")
                                    st.write(f"**Abstract:** {paper['abstract']}")
                                    if paper.get('full_text'):
                                        with st.expander("Show Full Text"):
                                            st.text(paper['full_text'])
                    
    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
