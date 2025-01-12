"""
Streamlit app for analyzing TdP risk of Ivabradine
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
import sys
import os
import json
import logging

# Add the parent directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
    dose = st.number_input("Dose (mg):", value=5.0, min_value=0.1, max_value=1000.0, step=0.1)
    analyze_button = st.button("Analyze")

# Main content
if analyze_button:
    try:
        with st.spinner("Analyzing drug..."):
            logger.info(f"Starting analysis for {drug_name}")
            
            # Get drug analysis
            analysis = st.session_state.analyzer.analyze_drug(drug_name, dose)
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
                    st.subheader("Drug Analysis")
                    
                    # CredibleMeds Box
                    with st.expander("CredibleMeds Status", expanded=True):
                        credible_meds = analysis.get('credible_meds', {})
                        if credible_meds:
                            st.write(f"**Category:** {credible_meds.get('category', 'Unknown')}")
                            st.write(f"**Risk Level:** {credible_meds.get('risk_level', 'Unknown')}")
                            if credible_meds.get('notes'):
                                st.write("**Notes:**")
                                st.write(credible_meds['notes'])
                        else:
                            st.write("No CredibleMeds data available")
                    
                    # hERG Binding Box
                    with st.expander("hERG Binding Prediction", expanded=True):
                        herg = analysis.get('herg', {})
                        if herg:
                            col1, col2 = st.columns(2)
                            with col1:
                                st.metric("IC50 Value", f"{herg.get('ic50', 'Unknown')} μM")
                            with col2:
                                st.metric("Safety Ratio", f"{herg.get('safety_ratio', 'Unknown')}x")
                                
                            if herg.get('notes'):
                                st.write("**Notes:**")
                                st.write(herg['notes'])
                        else:
                            st.write("No hERG binding data available")
                    
                    # Drug Concentration Box
                    with st.expander("Drug Concentrations", expanded=True):
                        conc = analysis.get('concentrations', {})
                        if conc:
                            col1, col2 = st.columns(2)
                            with col1:
                                st.metric("Therapeutic Concentration", f"{conc.get('therapeutic', 'Unknown')} μM")
                            with col2:
                                st.metric("Maximum Concentration", f"{conc.get('max', 'Unknown')} μM")
                        else:
                            st.write("No concentration data available")
                    
                    # Risk Assessment Box
                    with st.expander("Risk Assessment", expanded=True):
                        if analysis.get('crediblemeds_risk'):
                            risk_text = analysis.get('risk_category', 'Known Risk')
                            st.error(f"⚠️ {risk_text} (CredibleMeds)")
                            st.markdown("[View on CredibleMeds](https://crediblemeds.org)")
                        elif analysis.get('theoretical_binding'):
                            st.warning("⚠️ Potential hERG binding detected")
                        else:
                            st.success("✅ No significant hERG binding predicted")
                
                with tab2:
                    st.subheader("Literature Review")
                    literature = analysis.get('literature', {})
                    logger.info(f"Literature data: {literature}")
                    
                    if 'error' in literature:
                        st.error(literature['error'])
                    elif 'case_reports' in literature:
                        papers = literature['case_reports']
                        if isinstance(papers, list) and len(papers) > 0:
                            st.write(f"Found {len(papers)} relevant papers")
                            
                            # Display each paper in an expander
                            for i, paper in enumerate(papers, 1):
                                with st.expander(f"Case Report {i}: {paper.get('Title', 'Untitled Paper')}", expanded=False):
                                    col1, col2 = st.columns(2)
                                    
                                    with col1:
                                        st.markdown("**Patient Information**")
                                        st.write(f"Age: {paper.get('Age', 'N/A')}")
                                        st.write(f"Sex: {paper.get('Sex', 'N/A')}")
                                        
                                        if paper.get('Medical History'):
                                            st.markdown("**Medical History**")
                                            st.write(paper['Medical History'])
                                        
                                        if paper.get('Medication History'):
                                            st.markdown("**Medication History**")
                                            st.write(paper['Medication History'])
                                    
                                    with col2:
                                        st.markdown("**Medication Details**")
                                        if paper.get('Oral Dose (mg)'):
                                            st.write(f"Dose: {paper['Oral Dose (mg)']} mg")
                                        else:
                                            st.write("No dosing information available")
                                        
                                        st.markdown("**ECG Measurements**")
                                        if paper.get('Uncorrected QT (ms)'):
                                            st.write(f"QT: {paper['Uncorrected QT (ms)']} ms")
                                        if paper.get('QTc'):
                                            st.write(f"QTc: {paper['QTc']} ms")
                                        if paper.get('Heart Rate (bpm)'):
                                            st.write(f"Heart Rate: {paper['Heart Rate (bpm)']} bpm")
                                        if not any([paper.get('Uncorrected QT (ms)'), paper.get('QTc')]):
                                            st.write("No ECG measurements available")
                                        
                                        if paper.get('Blood Pressure (mmHg)'):
                                            st.write(f"Blood Pressure: {paper['Blood Pressure (mmHg)']} mmHg")
                                        
                                        if paper.get('Torsades de Pointes?') == 'Yes':
                                            st.error("⚠️ Torsades de Pointes reported")
                                    
                                    if paper.get('Course of Treatment'):
                                        st.markdown("**Treatment Course**")
                                        st.write(paper['Course of Treatment'])
                                    
                                    # Add source information
                                    st.markdown("---")
                                    source_text = []
                                    if paper.get('PMID'):
                                        source_text.append(f"PMID: {paper['PMID']}")
                                    if paper.get('DOI'):
                                        source_text.append(f"DOI: {paper['DOI']}")
                                    if source_text:
                                        st.caption(" | ".join(source_text))
                        else:
                            st.info("No case reports found matching the search criteria.")
                    else:
                        st.info("No literature data available.")

    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
