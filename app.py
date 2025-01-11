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
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Drug Information")
                        st.write(f"**PubChem CID:** {analysis.get('cid', 'Not found')}")
                        st.write(f"**Molecular Weight:** {analysis.get('molecular_weight', 'Not found')} g/mol")
                        st.write(f"**hERG IC50:** {analysis.get('herg_ic50', 'None')} μM")
                        st.write(f"**Source:** {analysis.get('source', 'Unknown')}")
                        
                        if 'theoretical_max' in analysis:
                            st.subheader("Concentration Analysis")
                            st.write(f"**Dose:** {analysis['dose']} mg")
                            st.write(f"**Theoretical Max:** {analysis['theoretical_max']:.2f} μM")
                            st.write(f"**Plasma Concentration:** {analysis['plasma_concentration']:.2f} μM")
                            if analysis.get('ratio_theoretical'):
                                st.write(f"**Theoretical/IC50 Ratio:** {analysis['ratio_theoretical']:.2f}")
                            if analysis.get('ratio_plasma'):
                                st.write(f"**Plasma/IC50 Ratio:** {analysis['ratio_plasma']:.2f}")
                    
                    with col2:
                        st.subheader("Risk Assessment")
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
                    
                    if 'error' in literature:
                        st.error(f"Error analyzing literature: {literature['error']}")
                    else:
                        # Show analysis results first
                        analysis_results = literature.get('analysis', [])
                        if analysis_results:
                            st.subheader("Case Report Analysis")
                            df = pd.DataFrame(analysis_results)
                            st.dataframe(df, use_container_width=True)
                            
                            # Show statistics
                            st.subheader("Summary Statistics")
                            col1, col2, col3 = st.columns(3)
                            
                            with col1:
                                tdp_cases = len([x for x in analysis_results if x.get('had_tdp') == 'Yes'])
                                st.metric("TdP Cases", tdp_cases)
                                
                            with col2:
                                avg_qtc = pd.DataFrame(analysis_results)['qtc'].mean()
                                if not pd.isna(avg_qtc):
                                    st.metric("Average QTc", f"{avg_qtc:.0f} ms")
                                    
                            with col3:
                                success_rate = len([x for x in analysis_results if x.get('treatment_successful') == 'Yes']) / len(analysis_results) * 100
                                st.metric("Treatment Success Rate", f"{success_rate:.0f}%")
                        
                        # Show papers by type
                        st.subheader("Literature Search Results")
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.subheader("Case Reports")
                            if literature.get('case_reports'):
                                df = pd.DataFrame(literature['case_reports'])
                                st.dataframe(df, use_container_width=True)
                            else:
                                st.write("No case reports found")
                        
                        with col2:
                            st.subheader("Cohort Studies")
                            if literature.get('cohort_studies'):
                                df = pd.DataFrame(literature['cohort_studies'])
                                st.dataframe(df, use_container_width=True)
                            else:
                                st.write("No cohort studies found")
                        
                        with col3:
                            st.subheader("Clinical Trials")
                            if literature.get('clinical_trials'):
                                df = pd.DataFrame(literature['clinical_trials'])
                                st.dataframe(df, use_container_width=True)
                            else:
                                st.write("No clinical trials found")
                    
    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
