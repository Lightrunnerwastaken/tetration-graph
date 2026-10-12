"""Read-only verification of a published research folder."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('folder',nargs='?',default=str(Path(__file__).absolute().parent))
a=p.parse_args();root=Path(a.folder).absolute()
manifest=json.loads((root/'MANIFEST_SHA256.json').read_text(encoding='utf-8'))
fail=[]
for r in manifest['files']:
    path=root/r['path']
    if not path.is_file():fail.append({'path':r['path'],'error':'missing'});continue
    data=path.read_bytes()
    if len(data)!=r['bytes'] or hashlib.sha256(data).hexdigest()!=r['sha256']:
        fail.append({'path':r['path'],'error':'hash or size mismatch'})
print(json.dumps({'root':str(root),'checked_payload_files':len(manifest['files']),
                  'payload_bytes':sum(r['bytes'] for r in manifest['files']),'failures':fail},indent=2))
raise SystemExit(bool(fail))
