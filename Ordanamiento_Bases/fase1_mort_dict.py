import pdfplumber

PDF_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Data dictionary for MCPS data file - mortality.pdf"

with pdfplumber.open(PDF_PATH) as pdf:
    total = len(pdf.pages)
    print(f"Total páginas: {total}")
    for pg_num, page in enumerate(pdf.pages, 1):
        text = page.extract_text()
        if text:
            print(f"\n{'='*60}\nPÁGINA {pg_num}\n{'='*60}\n{text}")
