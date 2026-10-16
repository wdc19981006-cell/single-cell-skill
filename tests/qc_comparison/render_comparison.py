"""Render historical/fixed PDFs with Poppler and build labelled comparison PNGs."""
import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader

PANELS = [('before_QC', 'before_QC'), ('after_QC', 'after_QC'),
          ('doublet_dimplot', 'doublet_umap'), ('doublet_vlnplot', 'doublet_vlnplot'),
          ('Singlet_dimplot', 'singlet_dimplot')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--fixed', type=Path)
    parser.add_argument('--dataset', action='append')
    args = parser.parse_args()
    poppler = shutil.which('pdftoppm')
    if not poppler:
        raise RuntimeError('Existing Poppler pdftoppm is required')
    args.output.mkdir(parents=True, exist_ok=True)
    font = ImageFont.truetype('C:/Windows/Fonts/arial.ttf', 24)
    validation = []
    for gse in args.dataset or ['GSE149614', 'GSE167297', 'GSE183904', 'GSE189903']:
        for original, skill in PANELS:
            paths = [('Original', args.reference / gse / (original + '.pdf')),
                     ('Old Skill', args.reference / gse / 'qc' / (skill + '.pdf'))]
            if args.fixed and gse == 'GSE167297':
                paths.append(('Fixed Skill', args.fixed / (skill + '.pdf')))
            width, height = 1600, 900
            canvas = Image.new('RGB', (width * len(paths), height + 80), 'white')
            draw = ImageDraw.Draw(canvas)
            for i, (label, pdf) in enumerate(paths):
                draw.text((i*width+20, 12), f'{gse} | {label} | {skill}', fill='black', font=font)
                if not pdf.exists():
                    draw.text((i*width+100, 400), 'Not produced by this version', fill='gray', font=font)
                    continue
                reader = PdfReader(pdf)
                assert len(reader.pages) == 1
                prefix = args.output / f'.render_{gse}_{i}_{skill}'
                subprocess.run([poppler, '-f', '1', '-singlefile', '-scale-to', '1600', '-png',
                                str(pdf), str(prefix)], check=True, capture_output=True)
                png = prefix.with_suffix('.png')
                with Image.open(png) as source:
                    source.load()
                    source.thumbnail((width, height))
                    canvas.paste(source, (i*width+(width-source.width)//2, 65+(height-source.height)//2))
                png.unlink()
                validation.append(dict(dataset=gse,version=label,plot=skill,pages=1,status='RENDER_PASS'))
            canvas.save(args.output / f'{gse}_{skill}.png', optimize=True)
            print('RENDERED', gse, skill, flush=True)
    with (args.output / 'render_validation.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=validation[0].keys())
        writer.writeheader()
        writer.writerows(validation)


if __name__ == '__main__':
    main()
