"""Open and render every final PDF/PNG pair; verify dimensions and 300 dpi."""
import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from PIL import Image
from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('qc', type=Path)
    parser.add_argument('renders', type=Path)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    args.renders.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in ('before_QC','after_QC','doublet_umap','doublet_vlnplot',
                 'singlet_dimplot','cell_cycle_phase','cell_cycle_score'):
        pdf = args.qc / (name+'.pdf')
        reader = PdfReader(pdf)
        assert len(reader.pages) == 1
        box = reader.pages[0].mediabox
        with Image.open(args.qc / (name+'.png')) as png:
            png.load()
            dpi = png.info['dpi']
            assert all(abs(x-300)<0.02 for x in dpi)
            inches = (float(box.width)/72, float(box.height)/72)
            assert all(abs(pixels-inch*300)<=2 for pixels,inch in zip(png.size,inches))
            rows.append(dict(plot=name,pdf_pages=1,width_inches=inches[0],height_inches=inches[1],
                             png_width=png.width,png_height=png.height,dpi=dpi[0],status='PASS'))
        subprocess.run([shutil.which('pdftoppm'),'-f','1','-singlefile','-scale-to','1600','-png',
                        str(pdf),str(args.renders/name)],check=True,capture_output=True)
    with args.report.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys())
        writer.writeheader();writer.writerows(rows)
    print('PASS: 7 PDFs rendered and 7 PNGs decoded; physical dimensions and 300 dpi verified')


if __name__ == '__main__':
    main()
