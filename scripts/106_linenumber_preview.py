"""A browser-viewable preview of the manuscript with visible line numbers -> docs/submission/line_number_preview.html
WHY THIS EXISTS. The submission .docx carries REAL Word line numbering (w:lnNumType in the section properties), which
Word, LibreOffice and Google Docs render in the margin and renumber as the author edits. That is the right thing for
submission and it cannot be seen without one of those applications. This file is a preview for everyone else: it
takes the same paragraphs out of the .docx and lays them out with numbered lines, so the numbering can be checked in
a browser.
WHAT IT IS NOT. The wrap points here are this script's, not Word's, so a given sentence will not necessarily fall on
the same numbered line in both. Nothing is submitted from this file and it is not the manuscript - it exists so that
"are the line numbers visible" can be answered by looking, on a machine with no Word.
Usage: python3 scripts/106_linenumber_preview.py
"""
import html
import textwrap
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs/submission/Periplasmic_surface_charge_manuscript.docx"
OUT = ROOT / "docs/submission/line_number_preview.html"
WRAP = 96          # characters per line, close to 6.5 inches of 12pt Times New Roman

doc = Document(SRC)
rows, n = [], 0
for p in doc.paragraphs:
    txt = p.text.rstrip()
    is_head = p.style.name.startswith("Heading")
    if not txt:
        n += 1
        rows.append((n, "", False))
        continue
    for ln in textwrap.wrap(txt, WRAP) or [""]:
        n += 1
        rows.append((n, ln, is_head))

body = "\n".join(
    f'<div class="l"><span class="n">{i}</span>'
    f'<span class="t{" h" if h else ""}">{html.escape(t) or "&nbsp;"}</span></div>'
    for i, t, h in rows)

OUT.write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Line number preview</title>
<style>
  :root {{ --ink:#111; --muted:#6b6b6b; --rule:#d8d8d4; --bg:#fcfcfb; }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{ --ink:#f0efec; --muted:#9a9a93; --rule:#3a3a38; --bg:#1a1a19; }}
  }}
  html,body {{ background:var(--bg); color:var(--ink); margin:0; }}
  body {{ font-family:"Times New Roman",Times,serif; font-size:16px; }}
  .page {{ max-width:60rem; margin:0 auto; padding:2rem 1rem 4rem; }}
  .note {{ font-family:system-ui,-apple-system,sans-serif; font-size:.82rem; color:var(--muted);
           border-left:3px solid var(--rule); padding:.5rem .9rem; margin:0 0 1.6rem; line-height:1.5; }}
  .l {{ display:flex; gap:1rem; align-items:baseline; }}
  .n {{ flex:0 0 3.2rem; text-align:right; font-family:system-ui,-apple-system,sans-serif;
        font-size:.72rem; color:var(--muted); user-select:none; border-right:1px solid var(--rule);
        padding-right:.5rem; }}
  .t {{ flex:1; line-height:2; white-space:pre-wrap; }}
  .t.h {{ font-weight:700; }}
  @media (max-width:640px) {{ .page {{ padding:1rem; }} .n {{ flex-basis:2.4rem; }} }}
</style></head><body><div class="page">
<p class="note"><strong>Preview only.</strong> The submission file is
<code>Periplasmic_surface_charge_manuscript.docx</code>, which carries real Word line numbering in its section
properties and renumbers itself as you edit. This page re-wraps the same text at a fixed width so the numbering can
be seen without Word; line breaks here will not match Word's exactly.</p>
{body}
</div></body></html>""")
print(f"wrote {OUT}  ({n} numbered lines)")
