import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import os
from ivablib.herg_analyzer import DrugAnalyzer

# Page config
st.set_page_config(
    page_title="TdP Risk Assessment",
    page_icon="❤️",
    layout="wide"
)

st.title("TdP Risk Assessment Tool")

# Load NCBI credentials from environment
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "your.email@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")

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
    # Simplified mock data
    risk_categories = {
        "ivabradine": "Known Risk of TdP",
        "amiodarone": "Known Risk of TdP",
        "sotalol": "Known Risk of TdP",
        "methadone": "Known Risk of TdP"
    }
    return risk_categories.get(drug_name, "Not Found in Database")

def create_risk_gauge(risk_category):
    """Create a Plotly gauge chart for risk visualization."""
    risk_levels = {
        "Known Risk of TdP": 1.0,
        "Possible Risk of TdP": 0.66,
        "Conditional Risk of TdP": 0.33,
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

if "hERG Channel Activity" in analysis_sections:
    st.markdown("## hERG Channel Activity", unsafe_allow_html=True)
    
    try:
        # Standard doses for ivabradine (mg)
        doses = [5.0, 7.5]  # Standard doses
        analysis = analyzer.analyze_drug(drug_name, doses)
        
        # Display hERG data
        st.markdown(f"""
        <div style='color: black;'>
        <p><b>Molecular Weight:</b> {analysis.molecular_weight:.2f} g/mol</p>
        <p><b>hERG IC50:</b> {analysis.herg_ic50 if analysis.herg_ic50 else 'Not found'} μM</p>
        <p><b>Source:</b> {analysis.herg_source if analysis.herg_source else 'Not available'}</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Display concentration data
        if analysis.concentrations:
            st.markdown("<div style='color: black;'><h3>Concentration Analysis</h3></div>", unsafe_allow_html=True)
            for i, conc in enumerate(analysis.concentrations):
                st.markdown(f"""
                <div style='color: black;'>
                <p><b>Dose {doses[i]} mg:</b></p>
                <ul>
                    <li>Theoretical Max: {conc.theoretical_max:.2f} μM</li>
                    <li>Plasma Concentration: {conc.plasma_concentration:.2f} μM</li>
                    {f'<li>IC50/Plasma Ratio: {conc.ratio_plasma:.2f}</li>' if conc.ratio_plasma else ''}
                </ul>
                </div>
                """, unsafe_allow_html=True)
        
        # Display citations
        st.markdown("<div style='color: black;'><h3>References</h3></div>", unsafe_allow_html=True)
        for citation in analysis.citations:
            st.markdown(f"<div style='color: black;'>{citation}</div>", unsafe_allow_html=True)
            
    except Exception as e:
        st.error(f"Error analyzing hERG activity: {str(e)}")
    
    st.markdown("## hERG Binding Analysis", unsafe_allow_html=True)
    st.markdown("<div style='color: black;'>Visualization of hERG channel binding characteristics will be added in future updates.</div>", unsafe_allow_html=True)

if "Literature Analysis" in analysis_sections:
    st.header("Literature Analysis")
    st.write("Literature analysis functionality will be added in future updates.")
