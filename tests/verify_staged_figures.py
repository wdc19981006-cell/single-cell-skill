"""Render actual stage PDFs with Poppler and verify PNG DPI and readable pages."""
import argparse
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path
from PIL import Image, ImageChops, ImageDraw


def run(project, output):
    project, output = Path(project).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    poppler = shutil.which("pdftoppm")
    pdfinfo = shutil.which("pdfinfo")
    if not poppler or not pdfinfo:
        raise RuntimeError("Poppler pdftoppm/pdfinfo required for actual PDF rendering")
    rows, previews = [], []
    files = sorted(p for p in project.rglob("*") if p.suffix.lower() in (".pdf", ".png")
        and not any(part in p.parts for part in (".workflow", ".runtime")) and not p.is_relative_to(output))
    for path in files:
        relative = path.relative_to(project).as_posix()
        if path.suffix.lower() == ".png":
            with Image.open(path) as im:
                im.load()
                dpi = im.info.get("dpi", (0, 0))
                if any(abs(value-300)>1 for value in dpi):
                    raise ValueError(f"Expected 300 DPI PNG: {relative}: {dpi}")
                rows.append(dict(file=relative,kind="PNG",pages=1,width=im.width,height=im.height,dpi=round(dpi[0]),status="PASS"))
            continue
        info = subprocess.run([pdfinfo,str(path)],check=True,capture_output=True,text=True,encoding="utf-8",errors="replace").stdout
        count = int(re.search(r"^Pages:\s+(\d+)",info,re.M).group(1))
        stem = output / (relative.replace("/", "__").removesuffix(".pdf"))
        subprocess.run([poppler,"-r","45","-png",str(path),str(stem)],check=True,capture_output=True)
        rendered = sorted(output.glob(stem.name+"-*.png"))
        if len(rendered)!=count:
            raise ValueError(f"Rendered page count mismatch: {relative}")
        for page in rendered:
            with Image.open(page) as im:
                rgb = im.convert("RGB")
                if ImageChops.difference(rgb,Image.new("RGB",rgb.size,"white")).getbbox() is None:
                    raise ValueError(f"Blank PDF page: {page}")
            previews.append((relative+" / "+page.stem.rsplit("-",1)[-1],page))
        rows.append(dict(file=relative,kind="PDF",pages=count,width="",height="",dpi="vector",status="PASS"))
    for page in range(0,len(previews),9):
        batch = previews[page:page+9]
        sheet = Image.new("RGB",(1200,1020),"#e3e6e9")
        draw = ImageDraw.Draw(sheet)
        for index,(label,path) in enumerate(batch):
            x,y=(index%3)*400,(index//3)*340
            with Image.open(path) as source:
                source.thumbnail((390,300))
                sheet.paste(source,(x+(400-source.width)//2,y+30))
            draw.text((x+6,y+5),label[-61:],fill="black")
        sheet.save(output/f"contact-{page//9+1}.jpg",quality=90)
    with (output/"figure_validation.csv").open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=("file","kind","pages","width","height","dpi","status"))
        writer.writeheader(); writer.writerows(rows)
    result=dict(status="PASS",files=len(rows),pdf_pages=len(previews),contact_sheets=(len(previews)+8)//9)
    (output/"test_results.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args(); run(args.project,args.output)
