"""Main application entry point."""
from ivablib import ui, pubmed, nlp, calculations, visualization
import os

def main():
    # Get user input
    drug_name = ui.get_drug_input()
    if not drug_name:
        ui.display_info("No drug name provided")
        return
        
    dose_mg = ui.get_dose_input()
    if not dose_mg:
        ui.display_info("No dose provided")
        return
    
    # Calculate concentrations
    ui.display_loading("Calculating drug concentrations")
    conc_results = calculations.calculate_drug_concentrations(drug_name, dose_mg)
    
    if not conc_results['theoretical_conc']:
        ui.display_error("Could not calculate concentrations. Missing molecular data.")
        return
    
    # Display concentration results
    ui.display_concentration_results(conc_results)
    ui.display_drug_info(conc_results)
    
    # Search and analyze literature
    papers = []
    qt_data = []
    processed_papers = []
    
    ui.display_loading("Searching medical literature")
    papers = pubmed.search_pubmed(drug_name)
    
    if not papers:
        ui.display_warning("No relevant papers found in literature search.")
        return
    
    ui.display_loading(f"Analyzing {len(papers)} papers for QT/HR data")
    for paper in papers:
        # Get full text if available
        text = paper.get('abstract', '') + ' ' + paper.get('title', '')
        full_text = pubmed.get_paper_full_text(paper)
        if full_text:
            text += ' ' + full_text
        
        # Extract values
        values = nlp.extract_qt_hr_values(text)
        if values['qt_values'] and values['hr_values']:
            for qt, _, _ in values['qt_values']:
                for hr, _, _ in values['hr_values']:
                    qt_data.append({'qt': qt, 'hr': hr})
        
        # Extract case details
        details = nlp.extract_case_details(text)
        paper.update(details)
        processed_papers.append(paper)
    
    # Display literature results
    ui.display_literature_results(
        papers_found=len(papers),
        qt_points=len([p for p in processed_papers if p.get('QT_values')]),
        hr_points=len([p for p in processed_papers if p.get('HR_values')])
    )
    
    # Create and save visualizations
    ui.display_loading("Generating visualizations")
    
    # Create output directory if it doesn't exist
    output_dir = "analysis_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save nomogram
    nomogram = visualization.create_nomogram(qt_data, drug_name)
    nomogram_path = os.path.join(output_dir, f"{drug_name}_nomogram.html")
    nomogram.write_html(nomogram_path)
    
    # Save analysis table
    table = visualization.create_analysis_table(processed_papers)
    table_path = os.path.join(output_dir, f"{drug_name}_analysis.html")
    table.write_html(table_path)
    
    ui.display_success(f"Analysis complete! Results saved to:")
    print(f"- Nomogram: {nomogram_path}")
    print(f"- Analysis Table: {table_path}")

if __name__ == "__main__":
    main()
