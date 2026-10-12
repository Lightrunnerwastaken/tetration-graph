"""Publish the prepared own research bundle to the chat output directory.

The explicitly requested D: research copy is made separately with native
PowerShell using sandbox escalation; this program writes only inside C: cwd.
"""
import hashlib,json,platform,shutil,sys
from pathlib import Path
import mpmath,numpy,matplotlib
HERE=Path(__file__).absolute().parent
TARGET=HERE.parents[1]/'outputs/hyper_halbrang'
origins={'halb_rang_hyperoperation_forschungsplan.md':r'D:\Users\janis\Downloads\halb_rang_hyperoperation_forschungsplan.md',
         'eta_values.csv':r'C:\Users\janis\Documents\Codex\2026-10-11\le\work\hyper_eta_rang\inputs\eta_values.csv',
         'previous_analysis.json':r'C:\Users\janis\Documents\Codex\2026-10-11\le\work\hyper_eta\analysis.json'}
provenance=[]
for name,origin in origins.items():
    data=(HERE/'inputs'/name).read_bytes();source=Path(origin)
    if source.read_bytes()!=data:raise ValueError(f'input differs from original: {name}')
    provenance.append({'file':name,'source':origin,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
(HERE/'inputs/provenance.json').write_text(json.dumps({'date':'2026-10-12','inputs':provenance,
    'python':sys.version,'platform':platform.platform(),'mpmath':mpmath.__version__,
    'numpy':numpy.__version__,'matplotlib':matplotlib.__version__,
    'upstream_snapshot_sha256':hashlib.sha256((HERE/'runtime/engine_snapshot.gp').read_bytes()).hexdigest(),
    'upstream_projects_modified':False,'shared_caches_used':False},indent=2),encoding='utf-8')
payload=[]
for path in sorted(HERE.rglob('*')):
    if not path.is_file() or '__pycache__' in path.parts or path.name.startswith('engine_with_hyper_'):continue
    if path.name=='MANIFEST_SHA256.json':continue
    rel=path.relative_to(HERE)
    data=path.read_bytes();dest=TARGET/rel;dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_bytes(data)
    payload.append({'path':rel.as_posix(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
manifest={'date':'2026-10-12','purpose':'Own half-rank research, not an established H_3.5',
          'files':payload,'excludes':['manifest itself','Python bytecode','regenerable combined GP scripts']}
(TARGET/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps({'target':str(TARGET),'payload_files':len(payload),'bytes':sum(r['bytes'] for r in payload)},indent=2))
