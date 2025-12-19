
Each PDF in this folder contains the source code as text from the `/source/` folder. Each *top-level* subfolder in `/source/` has its own PDF file. In the case of `IO`, each subfolder in `/source/IO/` has its own PDF file as well.

The source code is in MATLAB format and is intended to be easily searchable and copyable. This can be used in AI-based code assistants that can read PDF files to create a CryoGrid-specific chatbot. 

## Generating the PDFs

This repo includes a small Python script that regenerates the PDFs from the current `/source/` tree:

- Script: `pdf_code/generate_pdf_code.py`
- Dependency: `reportlab` (declared in `pdf_code/requirements.txt`)

Example (from repo root):

```bash
python pdf_code/generate_pdf_code.py --dry-run
python pdf_code/generate_pdf_code.py --overwrite
```

### Notes

- The script writes PDFs with embedded text (not raster images), so they remain searchable.
- Existing PDFs are left untouched unless you pass `--overwrite`.
