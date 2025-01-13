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
            
            # Get literature analysis using pubmed4125 module directly
            try:
                literature = analyze_literature(drug_name)
                if isinstance(literature, pd.DataFrame):
                    analysis['literature'] = literature.to_dict('records')
                else:
                    analysis['literature'] = literature
                logger.info(f"Literature analysis complete")
            except Exception as e:
                logger.error(f"Error in literature analysis: {e}")
                analysis['literature'] = {"error": str(e)}

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
                    
                    # Drug Concentrations Box
                    with st.expander("Drug Concentrations", expanded=True):
                        conc = analysis.get('concentrations', {})
                        if conc:
                            col1, col2 = st.columns(2)
                            with col1:
                                st.metric("Theoretical Max", f"{conc.get('theoretical_max', 'Unknown')} μM")
                                st.metric("Plasma Concentration", f"{conc.get('plasma_concentration', 'Unknown')} μM")
                            with col2:
                                st.metric("Theoretical/IC50 Ratio", conc.get('ratio_theoretical', 'Unknown'))
                                st.metric("Plasma/IC50 Ratio", conc.get('ratio_plasma', 'Unknown'))
                        else:
                            st.write("No concentration data available")
                
                with tab2:
                    st.subheader("Literature Review")
                    literature_data = analysis.get('literature', [])
                    
                    if isinstance(literature_data, list) and literature_data:
                        for paper in literature_data:
                            with st.expander(paper.get('Title', 'Untitled'), expanded=False):
                                st.write(f"**Authors:** {paper.get('Authors', 'Unknown')}")
                                st.write(f"**Journal:** {paper.get('Journal', 'Unknown')}")
                                st.write(f"**Year:** {paper.get('Year', 'Unknown')}")
                                st.write(f"**PMID:** {paper.get('PMID', 'Unknown')}")
                                
                                if paper.get('QT_Values'):
                                    st.markdown("**QT Values:**")
                                    st.write(paper['QT_Values'])
                                
                                if paper.get('TdP_Present'):
                                    st.markdown("**TdP Present:**")
                                    st.write(paper['TdP_Present'])
                                    
                                if paper.get('Medical_History'):
                                    st.markdown("**Medical History:**")
                                    st.write(paper['Medical_History'])
                                    
                                if paper.get('Treatment_Course'):
                                    st.markdown("**Treatment Course:**")
                                    st.write(paper['Treatment_Course'])
                    elif isinstance(literature_data, dict) and 'error' in literature_data:
                        st.error(f"Error in literature analysis: {literature_data['error']}")
                    else:
                        st.info("No literature data available.")

    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
