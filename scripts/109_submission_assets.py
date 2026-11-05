"""Export four journal figure files and the journal-hosted analysis-data supplement."""
from pathlib import Path
import hashlib
import shutil
import xml.etree.ElementTree as ET
import zipfile
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'results/broad2/figures'
OUT = ROOT / 'docs/submission/figures'
OUT.mkdir(parents=True, exist_ok=True)
for n, name in enumerate(['fig1_compartment_contrast', 'fig2_optimum_decomposition', 'fig3_per_reference_bh'], 1):
    for ext in ['png', 'svg']:
        shutil.copyfile(FIG / f'{name}.{ext}', OUT / f'Figure_{n}.{ext}')

parts = ['fig4a_censoring_sensitivity', 'fig4b_censoring_headroom']
ns = 'http://www.w3.org/2000/svg'
ET.register_namespace('', ns)
roots = [ET.parse(FIG / f'{p}.svg').getroot() for p in parts]
boxes = [list(map(float, r.attrib['viewBox'].split())) for r in roots]
width = max(b[2] for b in boxes)
gap = 30
height = sum(b[3] for b in boxes) + gap * 3
root = ET.Element(f'{{{ns}}}svg', {'width': f'{width}pt', 'height': f'{height}pt', 'viewBox': f'0 0 {width} {height}'})
y = gap
for i, (part, box) in enumerate(zip(roots, boxes)):
    # Scope IDs and references to prevent collisions between Matplotlib panels.
    raw = ET.tostring(part, encoding='unicode')
    import re
    ids = re.findall(r'\bid="([^"]+)"', raw)
    for key in sorted(ids, key=len, reverse=True):
        raw = raw.replace(f'id="{key}"', f'id="p{i}_{key}"').replace(f'#{key}', f'#p{i}_{key}')
    child = ET.fromstring(raw)
    child.set('x', '0'); child.set('y', str(y))
    label = ET.SubElement(root, f'{{{ns}}}text', {'x': '4', 'y': str(y - 7), 'font-family': 'sans-serif', 'font-size': '14', 'font-weight': 'bold'})
    label.text = chr(65 + i)
    root.append(child)
    y += box[3] + gap
ET.ElementTree(root).write(OUT / 'Figure_4.svg', encoding='utf-8', xml_declaration=True)

images = [Image.open(FIG / f'{p}.png').convert('RGB') for p in parts]
pad = 180
canvas = Image.new('RGB', (max(im.width for im in images), sum(im.height for im in images) + 3 * pad), 'white')
draw = ImageDraw.Draw(canvas)
font = ImageFont.truetype(font_manager.findfont('DejaVu Sans'), 88)
y = pad
for i, im in enumerate(images):
    draw.text((20, y - 120), chr(65 + i), fill='black', font=font)
    canvas.paste(im, (0, y))
    y += im.height + pad
canvas.save(OUT / 'Figure_4.png', dpi=(600, 600))

# Original raw downloads remain external; this archive holds the curated inputs
# and stored numerical outputs needed to reproduce figures and inspect claims.
paths = []
for folder in ['data/curated', 'results/broad2', 'results/phylo_broad2']:
    paths.extend(p for p in (ROOT / folder).rglob('*') if p.is_file() and '/figures/' not in p.as_posix())
paths.extend(ROOT / p for p in ['environment.yml', 'requirements.txt', 'PROVENANCE.md', 'DATA_AVAILABILITY.md', 'DATA_LICENSE.md', 'LICENSE', 'docs/DATA_DICTIONARY.md', 'docs/NUMBERS.md', 'docs/RELEASE.md'])
paths.extend((ROOT / 'results').glob('*classes.csv'))
paths = sorted(set(paths))
archive = ROOT / 'docs/submission/Supplemental_Data_S1.zip'
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    sums = []
    for p in paths:
        name = p.relative_to(ROOT).as_posix()
        z.write(p, name)
        sums.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {name}')
    z.writestr('SUPPLEMENT_CHECKSUMS.sha256', '\n'.join(sums) + '\n')
    z.writestr('README.txt', 'Supplemental Data S1\nCurated inputs and stored analysis outputs for the extracytoplasmic charge study.\nSee DATA_AVAILABILITY.md for excluded third-party raw inputs, retrieval instructions and manual curation limits.\nAnalysis code and reproducibility commands: https://github.com/ajayvellanki0720-boop/acidophile-ph-limits-submission\nThis file is intended for direct upload as journal-hosted supplemental data, without an external DOI.\n')
assert archive.stat().st_size < 30_000_000, 'Split supplement before upload: exceeds 30 MB.'
print(f'Created four PNG/SVG figure pairs and {archive.stat().st_size / 1e6:.2f} MB supplement')
