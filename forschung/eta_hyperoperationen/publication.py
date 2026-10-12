"""Assemble the self-contained local research deliverable; no source project edits."""
import hashlib,json,platform,shutil,sys
from pathlib import Path
from common import HERE
target=Path(r'C:\Users\janis\Documents\Codex\2026-10-11\le\outputs\hyper_eta')
target.mkdir(parents=True,exist_ok=True)
allowed={'.md','.json','.csv','.png','.svg','.py','.gp','.log'}
archived={'probe_hex.py','probe_hex_2_25.json','hyper_real.gp'}
for f in sorted(HERE.iterdir()):
    if not f.is_file() or f.suffix not in allowed:continue
    dest=(target/'archiv'/f.name) if f.name in archived or f.name.startswith('engine_with_hyper') else target/f.name
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,dest)
for f in (HERE/'_backend/fatou_backend').glob('*.py'):
    dest=target/'_backend/fatou_backend'/f.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,dest)
(target/'archiv/README.md').write_text('Historische Arbeitsstände. Die GP-Kombinationen sichern die während der Untersuchung entstandenen Codeversionen. Der erste probe_hex-Hilfsprototyp hatte nach seiner gültigen Initialisierungsabfrage einen GP-Syntaxfehler; er wird nur als unveränderter Arbeitsstand gesichert. Empfohlene reproduzierbare Programme und Ergebnisse stehen im Hauptordner.\n',encoding='utf-8')
files=[]
for f in sorted(target.rglob('*')):
    if f.is_file() and f.name!='manifest.json':
        files.append({'path':str(f.relative_to(target)).replace('\\','/'),'bytes':f.stat().st_size,
                      'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
manifest={'date':'2026-10-12','python':platform.python_version(),'files':files,
          'destination_research_folder':r'D:\Users\janis\Documents\Matura\tetration-graph\forschung\eta_hyperoperationen',
          'proof_project_modified':False,'upstream_sources_and_shared_caches_modified':False}
(target/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'files':len(files),'bytes':sum(z['bytes'] for z in files),'target':str(target)},indent=2))
