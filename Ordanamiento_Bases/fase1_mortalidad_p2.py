import pdfplumber

PDF_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Listing - ICD-10 codes contributing to each derived MCPS mortality endpoint.pdf"

with pdfplumber.open(PDF_PATH) as pdf:
    # Páginas 2 y 3 tienen los códigos de cáncer
    for pg in [1, 2, 3]:
        text = pdf.pages[pg].extract_text()
        print(f"\n{'='*60}\nPÁGINA {pg+1}\n{'='*60}\n{text}")
