import pdfplumber
import re

PDF_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Data dictionary for MCPS data file - baseline.pdf"

TERMS = [
    "cancer", "tumour", "tumor", "malignancy", "neoplasm",
    "breast", "prostate", "lung", "gastrointestinal",
    "leukemia", "leukaemia", "lymphoma", "skin", "cervix",
    "uterus", "ovary", "colon", "rectum", "rectal",
    "carcinoma", "sarcoma", "melanoma", "oncolog"
]

pattern = re.compile("|".join(TERMS), re.IGNORECASE)

results = []
with pdfplumber.open(PDF_PATH) as pdf:
    total_pages = len(pdf.pages)
    print(f"Total páginas: {total_pages}")
    for page_num, page in enumerate(pdf.pages, 1):
        text = page.extract_text()
        if text and pattern.search(text):
            lines = text.split("\n")
            for i, line in enumerate(lines):
                if pattern.search(line):
                    context_start = max(0, i-1)
                    context_end = min(len(lines), i+3)
                    context = "\n".join(lines[context_start:context_end])
                    results.append({
                        "page": page_num,
                        "line": line.strip(),
                        "context": context.strip()
                    })

print(f"\nTotal coincidencias en líneas: {len(results)}\n")
for r in results:
    print(f"--- Pág {r['page']} ---")
    print(r['context'])
    print()
