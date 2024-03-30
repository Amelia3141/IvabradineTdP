import os
import PyPDF2
import matplotlib.pyplot as plt
import spacy
import re  # Regular expression library

# Load the spaCy English model with NER
nlp = spacy.load("en_core_web_sm")

def extract_text_from_pdf(pdf_path):
    with open(pdf_path, 'rb') as file:
        reader = PyPDF2.PdfReader(file)
        text = ""
        for page_num in range(len(reader.pages)):
            text += reader.pages[page_num].extract_text()
        return text

def extract_qt_hr_values(pdf_text):
    qt_hr_pairs = []
    doc = nlp(pdf_text)

    for ent in doc.ents:
        if ent.label_ == 'CARDINAL':
            value = ent.text
            if 'ms' in ent.sent.text:
                qt_value_pattern = re.compile(r'QT.?:?\s*(\d{3,4})\s*ms', re.IGNORECASE)
                qt_match = qt_value_pattern.search(ent.sent.text)
                if qt_match:
                    qt_value = int(qt_match.group(1))
                    for token in ent.sent:
                        if token.text.lower() == 'hr' or token.text.lower() == 'heart rate':
                            for child in token.children:
                                if child.ent_type_ == 'CARDINAL':
                                    hr_value_pattern = re.compile(r'(\d{2,3})\s*bpm', re.IGNORECASE)
                                    hr_match = hr_value_pattern.search(child.sent.text)
                                    if hr_match:
                                        hr_value = int(hr_match.group(1))
                                        qt_hr_pairs.append((qt_value, hr_value))
                                        break

    return qt_hr_pairs

def plot_qt_nomogram(qt_hr_pairs, case_titles):
    # Placeholder for the body of your plotting function
    # Insert your actual plotting code here.
    pass  # Remove this line and replace it with your plotting code

# Usage example
folder_path = '/Users/ameliag/ivabcode'
pdf_files = [file for file in os.listdir(folder_path) if file.endswith('.pdf')]

qt_hr_pairs = []
case_titles = []

for pdf_file in pdf_files:
    pdf_path = os.path.join(folder_path, pdf_file)
    if os.path.isfile(pdf_path):
        pdf_text = extract_text_from_pdf(pdf_path)
        qt_hr_values = extract_qt_hr_values(pdf_text)
        if qt_hr_values:
            qt_hr_pairs.extend(qt_hr_values)  # Use extend to flatten the list
            case_titles.append(pdf_file)
    else:
        print(f"Error: {pdf_path} is not a valid file.")

if qt_hr_pairs and case_titles:
    plot_qt_nomogram(qt_hr_pairs, case_titles)
else:
    print("No QT or heart rate values found in the PDF files.")
