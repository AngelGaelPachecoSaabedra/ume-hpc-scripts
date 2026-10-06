import pdfplumber

PDF_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Listing - ICD-10 codes contributing to each derived MCPS mortality endpoint.pdf"

with pdfplumber.open(PDF_PATH) as pdf:
    total_pages = len(pdf.pages)
    print(f"Total páginas: {total_pages}\n")
    for page_num, page in enumerate(pdf.pages, 1):
        text = page.extract_text()
        if text:
            print(f"\n{'='*60}")
            print(f"PÁGINA {page_num}")
            print('='*60)
            print(text)
