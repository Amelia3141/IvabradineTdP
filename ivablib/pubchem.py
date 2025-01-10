import requests
import json
import time

def get_herg_activity(drug_name):
    """Get hERG activity information from PubChem."""
    try:
        # Search for the compound
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
        response = requests.get(search_url)
        if response.status_code != 200:
            return "Unable to find drug in PubChem database"
            
        data = response.json()
        if "IdentifierList" not in data:
            return "No PubChem data available for this drug"
            
        cid = data["IdentifierList"]["CID"][0]
        
        # Add delay to respect PubChem's rate limits
        time.sleep(0.5)
        
        # Get bioassay data
        assay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
        response = requests.get(assay_url)
        if response.status_code != 200:
            return "Unable to retrieve bioassay data"
            
        data = response.json()
        
        # Look for hERG-related assays
        if "AssaySummaries" in data:
            herg_assays = []
            for assay in data["AssaySummaries"]["AssaySummary"]:
                assay_name = assay["Name"].lower()
                if any(term in assay_name for term in ["herg", "kcnh2", "potassium channel", "qt"]):
                    result = "Active" if assay.get("ActivityOutcomeMethod", "").lower() == "active" else "Tested"
                    herg_assays.append(f"{result} in {assay['Name']} (AID: {assay['AID']})")
            
            if herg_assays:
                return "hERG Activity Found:\n" + "\n".join(herg_assays)
        
        return "No specific hERG activity data found in PubChem"
        
    except Exception as e:
        print(f"Error in PubChem search: {e}")
        return "Error retrieving hERG activity data"
