"""Data visualization components."""
import plotly.graph_objects as go
from typing import Dict, List
import pandas as pd

def create_nomogram(qt_data: List[Dict[str, float]], drug_name: str) -> go.Figure:
    """Create QT nomogram plot."""
    # Reference lines data
    ref_lines = {
        'upper': [(0, 400), (60, 440), (120, 500)],
        'lower': [(0, 300), (60, 340), (120, 400)]
    }
    
    # Create figure
    fig = go.Figure()
    
    # Add reference lines
    for line_type, points in ref_lines.items():
        x, y = zip(*points)
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode='lines',
            name='At Risk' if line_type == 'upper' else 'Normal',
            line=dict(color='red' if line_type == 'upper' else 'green')
        ))
    
    # Add data points
    if qt_data:
        hr_values = [point['hr'] for point in qt_data]
        qt_values = [point['qt'] for point in qt_data]
        fig.add_trace(go.Scatter(
            x=hr_values,
            y=qt_values,
            mode='markers',
            name=drug_name,
            marker=dict(size=10)
        ))
    
    # Update layout
    fig.update_layout(
        title=f"QT Nomogram - {drug_name}",
        xaxis_title="Heart Rate (bpm)",
        yaxis_title="QT Interval (ms)",
        showlegend=True
    )
    
    return fig

def create_analysis_table(papers: List[Dict[str, str]]) -> go.Figure:
    """Create analysis results table."""
    # Convert papers to DataFrame
    df = pd.DataFrame(papers)
    
    # Select relevant columns
    display_cols = ['Title', 'Authors', 'Year', 'QT_values', 'HR_values', 'Medical_History']
    df = df[display_cols]
    
    # Create table
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=list(df.columns),
            fill_color='paleturquoise',
            align='left'
        ),
        cells=dict(
            values=[df[col] for col in df.columns],
            fill_color='lavender',
            align='left'
        )
    )])
    
    fig.update_layout(
        title="Literature Analysis Results",
        height=400
    )
    
    return fig
