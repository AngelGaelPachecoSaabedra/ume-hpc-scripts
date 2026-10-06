import pdfplumber

PDF_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Data dictionary for MCPS data file - baseline.pdf"

with pdfplumber.open(PDF_PATH) as pdf:
    for page_num, page in enumerate(pdf.pages, 1):
        text = page.extract_text()
        if text:
            print(f"\n{'='*60}")
            print(f"PÁGINA {page_num}")
            print('='*60)
            print(text)
