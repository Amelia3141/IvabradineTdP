"""
Streamlit app for analyzing TdP risk of Ivabradine
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import json
import logging
from typing import Dict, List, Optional, Tuple
from ivablib import DrugAnalyzer
from ivablib.case_report_analyzer import CaseReportAnalyzer

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
    
if 'case_analyzer' not in st.session_state:
    st.session_state.case_analyzer = CaseReportAnalyzer()
    
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
                        st.header("Drug Information")
                        
                        # Display hERG IC50 with source
                        herg_ic50 = analysis.get('herg_ic50')
                        herg_source = analysis.get('herg_source', 'No data available')
                        
                        # Check if hERG IC50 is a number before formatting
                        if isinstance(herg_ic50, (int, float)):
                            st.metric("hERG IC50", f"{herg_ic50:.2f} μM")
                        else:
                            st.metric("hERG IC50", str(herg_ic50))
                        
                        st.caption(f"Source: {herg_source}")
                            
                        # Add more drug information here if available
                        if analysis.get('drug_info'):
                            st.write("**Additional Information:**")
                            for key, value in analysis['drug_info'].items():
                                if key not in ['herg_ic50', 'herg_source']:
                                    st.write(f"**{key.replace('_', ' ').title()}:** {value}")
                    
                    with col2:
                        st.subheader("Risk Assessment")
                        if analysis.get('crediblemeds_risk'):
                            risk_text = analysis.get('risk_category', 'Known Risk of TdP')
                            st.error(f"⚠️ {risk_text.title()} (CredibleMeds)")
                            st.markdown("[View on CredibleMeds](https://crediblemeds.org)")
                        elif analysis.get('herg_ic50') and analysis['herg_ic50'] < 10:
                            st.warning("⚠️ Potential Risk: hERG IC50 < 10 μM")
                        else:
                            st.info("ℹ️ No known TdP risk data available")
                
                with tab2:
                    st.subheader("Literature Review")
                    literature = analysis.get('literature', pd.DataFrame())
                    if literature.empty:
                        st.info("No relevant literature found.")
                    else:
                        # Process case reports
                        case_reports = st.session_state.case_analyzer.analyze_papers(
                            literature.to_dict('records'),
                            drug_name
                        )
                        
                        if not case_reports.empty:
                            # Show summary table first
                            st.write("Summary of Case Reports:")
                            
                            # Create a display dataframe with all needed columns
                            display_df = case_reports[['title', 'journal', 'year', 'authors', 'age', 'sex', 'qtc', 'had_tdp']].copy()
                            display_df.columns = ['Title', 'Journal', 'Year', 'Authors', 'Age', 'Sex', 'QTc', 'TdP']
                            
                            # Truncate long titles and authors for better display
                            display_df['Title'] = display_df['Title'].str.slice(0, 100) + '...'
                            display_df['Authors'] = display_df['Authors'].str.slice(0, 50) + '...'
                            
                            # Add index starting from 1
                            display_df.index = range(1, len(display_df) + 1)
                            
                            # Display with full width and column config
                            st.dataframe(
                                display_df,
                                use_container_width=True,
                                column_config={
                                    "Title": st.column_config.TextColumn(width="large"),
                                    "Journal": st.column_config.TextColumn(width="medium"),
                                    "Year": st.column_config.NumberColumn(width="small"),
                                    "Authors": st.column_config.TextColumn(width="medium"),
                                    "Age": st.column_config.NumberColumn(width="small"),
                                    "Sex": st.column_config.TextColumn(width="small"),
                                    "QTc": st.column_config.NumberColumn(width="small"),
                                    "TdP": st.column_config.TextColumn(width="small")
                                }
                            )
                        
                        # Show detailed expandable list
                        st.write("\nDetailed Paper Information:")
                        for idx, paper in case_reports.iterrows():
                            with st.expander(f"{idx+1}. {paper['title']}"):
                                st.write(f"**Authors:** {paper['authors']}")
                                st.write(f"**Journal:** {paper['journal']} ({paper['year']})")
                                if paper['age'] or paper['sex']:
                                    st.write(f"**Patient:** {paper['age']} years old, {paper['sex']}")
                                if paper['qtc']:
                                    st.write(f"**QTc:** {paper['qtc']} ms")
                                if paper['had_tdp']:
                                    st.write(f"**TdP:** {paper['had_tdp']}")
                                if paper['treatment_successful']:
                                    st.write(f"**Treatment Outcome:** {paper['treatment_successful']}")
                                if paper['drug_combinations']:
                                    st.write(f"**Drug Combinations:** {paper['drug_combinations']}")
                                st.write(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/)")
    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
else:
    st.info("Click 'Analyze' to start the analysis")
