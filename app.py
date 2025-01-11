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
    st.session_state.drug_name = "ivabradine"
if 'dose' not in st.session_state:
    st.session_state.dose = 5.0  # Default dose in mg

# Sidebar inputs
with st.sidebar:
    st.header("Analysis Parameters")
    
    # Drug selection
    drug_name = st.text_input("Drug Name", value=st.session_state.drug_name)
    
    # Single dose input
    dose = st.number_input("Dose (mg)", 
                          value=st.session_state.dose,
                          min_value=0.0,
                          max_value=1000.0,
                          step=0.5)
    
    # Update session state
    st.session_state.drug_name = drug_name
    st.session_state.dose = dose
    
    analyze_button = st.button("Analyze")

def search_literature(drug_name: str) -> pd.DataFrame:
    """Search PubChem literature for drug information"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/xrefs/PubMedID/JSON"
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        if 'InformationList' in data and 'Information' in data['InformationList']:
            pmids = data['InformationList']['Information'][0].get('PubMedID', [])
            
            # Get details for each PubMed ID
            papers = []
            for pmid in pmids[:10]:  # Limit to first 10 papers
                try:
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
                    response = requests.get(url)
                    response.raise_for_status()
                    paper_data = response.json()
                    
                    if 'result' in paper_data and str(pmid) in paper_data['result']:
                        paper = paper_data['result'][str(pmid)]
                        papers.append({
                            'Title': paper.get('title', ''),
                            'Authors': ', '.join(paper.get('authors', [])),
                            'Journal': paper.get('source', ''),
                            'Year': paper.get('pubdate', '').split()[0],
                            'PMID': pmid
                        })
                    
                    time.sleep(0.1)  # Rate limiting
                    
                except Exception as e:
                    st.warning(f"Error fetching paper {pmid}: {str(e)}")
                    continue
            
            return pd.DataFrame(papers)
    except Exception as e:
        st.error(f"Error searching literature: {str(e)}")
        return pd.DataFrame()

# Main content
if analyze_button:
    try:
        with st.spinner("Analyzing drug..."):
            logger.info(f"Starting analysis for {drug_name}")
            analysis = st.session_state.analyzer.analyze_drug(drug_name, [dose])
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
                        st.write(f"**PubChem CID:** {analysis.get('cid', 'Not found')}")
                        st.write(f"**Molecular Weight:** {analysis.get('molecular_weight', 'Not found')} g/mol")
                        st.write(f"**hERG IC50:** {analysis.get('herg_ic50', 'None')} μM")
                        st.write(f"**Source:** {analysis.get('source', 'Unknown')}")
                    
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
                    
                    # Concentration analysis
                    st.subheader("Concentration Analysis")
                    
                    # Convert concentrations to DataFrame
                    df = pd.DataFrame(analysis['concentrations'])
                    df.columns = [
                        'Dose (mg)',
                        'Theoretical Max (μM)',
                        'Plasma Concentration (μM)',
                        'Theoretical/IC50 Ratio',
                        'Plasma/IC50 Ratio'
                    ]
                    
                    # Format numbers
                    for col in df.columns[1:]:
                        df[col] = df[col].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else "N/A")
                    
                    st.dataframe(df)
                    
                    # Visualization
                    st.subheader("Concentration vs IC50 Visualization")
                    
                    # Create traces for theoretical and plasma concentrations
                    fig = go.Figure()
                    
                    # Add IC50 line if available
                    if analysis['herg_ic50']:
                        fig.add_hline(y=analysis['herg_ic50'], 
                                    line_dash="dash", 
                                    line_color="red",
                                    annotation_text="hERG IC50")
                    
                    # Add concentration points
                    fig.add_trace(go.Scatter(
                        x=[dose],
                        y=[float(df['Theoretical Max (μM)'].iloc[0]) if df['Theoretical Max (μM)'].iloc[0] != "N/A" else None],
                        name="Theoretical Max",
                        mode="markers",
                        marker=dict(color='blue', size=10)
                    ))
                    
                    fig.add_trace(go.Scatter(
                        x=[dose],
                        y=[float(df['Plasma Concentration (μM)'].iloc[0]) if df['Plasma Concentration (μM)'].iloc[0] != "N/A" else None],
                        name="Plasma Concentration",
                        mode="markers",
                        marker=dict(color='green', size=10)
                    ))
                    
                    fig.update_layout(
                        title="Drug Concentrations vs hERG IC50",
                        xaxis_title="Dose (mg)",
                        yaxis_title="Concentration (μM)",
                        showlegend=True
                    )
                    
                    st.plotly_chart(fig)
                    
                    # Citations
                    st.subheader("References")
                    for citation in analysis['citations']:
                        st.markdown(f"- {citation}")
                
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
else:
    st.info("Click 'Analyze' to start the analysis")
