"""Create an immutable user-facing output package and SHA256 manifest."""
import json,hashlib,shutil
from pathlib import Path
HERE=Path(__file__).absolute().parent
DEST=HERE.parent.parent/'outputs'/'hyper_halbrang_analytisch'
def run():
    files=[]
    for p in HERE.rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts:continue
        rel=p.relative_to(HERE)
        if rel.parts[0] not in {'inputs','data','figures'} and p.suffix not in {'.py','.md'}:continue
        files.append((p,rel))
    if DEST.exists():raise FileExistsError(f'Publication folder already exists: {DEST}')
    DEST.mkdir(parents=True)
    records=[]
    for p,rel in sorted(files,key=lambda item:str(item[1])):
        out=DEST/rel;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
        records.append({'path':rel.as_posix(),'bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()})
    manifest={'date':'2026-10-12','scope':'analytic quadratic half-rank jet; no complete H3.5 or eta3.5',
        'payload_files':len(records),'payload_bytes':sum(r['bytes'] for r in records),'files':records}
    (DEST/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'destination':str(DEST),'payload_files':len(records),'payload_bytes':manifest['payload_bytes']}))
if __name__=='__main__':run()
