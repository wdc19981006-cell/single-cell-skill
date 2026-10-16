"""Render clearly labelled historical/control/current UMAP and PCA previews."""
import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('latest',type=Path)
    parser.add_argument('control',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',24)
    small=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
    validation=[]
    for stage,original,stem in [('pre_filter','doublet_dimplot','doublet_umap'),
                                ('post_filter','Singlet_dimplot','singlet_dimplot')]:
        panels=[('Original code (historical)',args.source/(original+'.pdf')),
                ('Legacy PC rule CONTROL (recreated)',args.control/'qc'/(stem+'.pdf')),
                ('Current main (independent run)',args.latest/'qc'/(stem+'.pdf'))]
        w,h=1800,1000
        canvas=Image.new('RGB',(w*3,h+130),'white'); draw=ImageDraw.Draw(canvas)
        for index,(label,pdf) in enumerate(panels):
            reader=PdfReader(pdf);assert len(reader.pages)==1
            draw.text((index*w+20,15),f'{stage} | {label}',fill='black',font=font)
            text=('Historical coordinates unavailable; genes/regression differ' if index==0 else
                  'Same raw and cell-QC input; pMT regression unchanged')
            draw.text((index*w+20,50),text,fill='#555555',font=small)
            if stage=='post_filter':
                draw.text((index*w+20,80),'Actual retained sets differ; see matched-cell comparison separately',fill='#555555',font=small)
            prefix=args.output/f'.render_{stage}_{index}'
            subprocess.run([shutil.which('pdftoppm'),'-f','1','-singlefile','-scale-to','2000','-png',
                            str(pdf),str(prefix)],check=True,capture_output=True)
            png=prefix.with_suffix('.png')
            with Image.open(png) as image:
                image.load();image.thumbnail((w,h))
                canvas.paste(image,(index*w+(w-image.width)//2,130+(h-image.height)//2))
            png.unlink()
            validation.append(dict(stage=stage,version=label,render='PASS'))
        canvas.save(args.output/f'{stage}_three_versions.png',optimize=True)
    elbows=sorted((args.latest/'diagnostics').glob('elbow_GSM*.png'))
    assert len(elbows)==20
    montage=Image.new('RGB',(2400,2500),'white')
    for index,path in enumerate(elbows):
        with Image.open(path) as image:
            image.load();image.thumbnail((600,500))
            montage.paste(image,((index%4)*600,(index//4)*500))
    montage.save(args.output/'sample_elbows_overview.png',optimize=True)
    with (args.output/'render_validation.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=validation[0].keys());writer.writeheader();writer.writerows(validation)
    print('PASS historical/control/current PDFs rendered and 20-sample elbow overview created')


if __name__=='__main__':
    main()
