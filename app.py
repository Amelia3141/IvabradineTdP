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
from typing import Dict, List, Optional, Tuple

from ivablib.herg_analyzer import DrugAnalyzer

# Page config
st.set_page_config(
    page_title="Ivabradine TdP Risk Analyzer",
    page_icon="💊",
    layout="wide"
)

# Title and description
st.title("Ivabradine TdP Risk Analyzer")
st.markdown("""
This app analyzes the Torsades de Pointes (TdP) risk of Ivabradine by:
1. Retrieving hERG IC50 values from PubChem
2. Calculating theoretical and plasma concentrations
3. Estimating TdP risk based on concentration ratios
""")

# Initialize session state
if 'drug_name' not in st.session_state:
    st.session_state.drug_name = "ivabradine"
if 'doses' not in st.session_state:
    st.session_state.doses = [5.0, 7.5, 10.0]  # Default doses in mg

# Sidebar inputs
with st.sidebar:
    st.header("Analysis Parameters")
    
    # Drug selection
    drug_name = st.text_input("Drug Name", value=st.session_state.drug_name)
    
    # Dose inputs
    st.subheader("Doses (mg)")
    doses = []
    for i, default_dose in enumerate(st.session_state.doses):
        dose = st.number_input(f"Dose {i+1}", 
                             value=default_dose,
                             min_value=0.0,
                             max_value=1000.0,
                             step=0.5)
        doses.append(dose)
    
    # Update session state
    st.session_state.drug_name = drug_name
    st.session_state.doses = doses
    
    analyze_button = st.button("Analyze")

# Main content
if analyze_button:
    try:
        with st.spinner("Analyzing drug..."):
            # Initialize analyzer with NCBI credentials
            analyzer = DrugAnalyzer(
                email=st.secrets["NCBI_EMAIL"],
                api_key=st.secrets["NCBI_API_KEY"]
            )
            
            # Analyze drug
            analysis = analyzer.analyze_drug(drug_name, doses)
            
            if "error" in analysis:
                st.error(analysis["error"])
            else:
                # Display results
                col1, col2 = st.columns(2)
                
                with col1:
                    st.subheader("Drug Information")
                    st.write(f"**Name:** {analysis['name']}")
                    st.write(f"**PubChem CID:** {analysis['cid']}")
                    st.write(f"**Molecular Weight:** {analysis['molecular_weight']} g/mol")
                    st.write(f"**hERG IC50:** {analysis['herg_ic50']} μM")
                    st.write(f"**Source:** {analysis['herg_source']}")
                
                with col2:
                    st.subheader("Risk Assessment")
                    if analysis['theoretical_binding']:
                        st.warning("⚠️ Potential hERG binding detected")
                    else:
                        st.success("✅ No significant hERG binding predicted")
                
                # Concentration table
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
                
                # Add concentration traces
                fig.add_trace(go.Scatter(
                    x=df['Dose (mg)'],
                    y=[float(x) if x != "N/A" else None for x in df['Theoretical Max (μM)']],
                    name="Theoretical Max",
                    line=dict(color='blue')
                ))
                
                fig.add_trace(go.Scatter(
                    x=df['Dose (mg)'],
                    y=[float(x) if x != "N/A" else None for x in df['Plasma Concentration (μM)']],
                    name="Plasma Concentration",
                    line=dict(color='green')
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
                    
    except Exception as e:
        st.error(f"Error analyzing drug: {str(e)}")
else:
    st.info("Click 'Analyze' to start the analysis")
