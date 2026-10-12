"""Publish this completed round and verify all copies; never remove files."""
import argparse
import ast
import hashlib
import json
import re
import shutil
from pathlib import Path

HERE=Path(__file__).absolute().parent
OUTPUT=Path(r"C:\Users\janis\Documents\Codex\2026-10-11\le\outputs\hyper_eta_rang")
RESEARCH=Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\eta_hyperoperationen\rang_gamma")
p=argparse.ArgumentParser()
p.add_argument("--research",action="store_true")
a=p.parse_args()

for py in HERE.rglob("*.py"):
    ast.parse(py.read_text(encoding="utf-8-sig"),filename=str(py))
for name in ["README.md","BERICHT.md","MATHEMATIK.md"]:
    prose=re.sub(r"\\\[.*?\\\]","",(HERE/name).read_text(encoding="utf-8"),flags=re.S)
    for match in re.finditer(r"\]\(([^)]+)\)",prose):
        ref=match.group(1)
        if "://" not in ref and not (HERE/ref).is_file():
            raise ValueError(f"broken relative link in {name}: {ref}")
files=sorted(f for f in HERE.rglob("*") if f.is_file() and f.name!="manifest.json" and "__pycache__" not in f.parts)
entries=[{"path":f.relative_to(HERE).as_posix(),"bytes":f.stat().st_size,
          "sha256":hashlib.sha256(f.read_bytes()).hexdigest()} for f in files]
manifest={"date":"2026-10-12","outputs":str(OUTPUT),"research_copy":str(RESEARCH),
          "file_count_excluding_manifest":len(entries),"entries":entries}
(HERE/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
files.append(HERE/"manifest.json")
targets=[RESEARCH] if a.research else [OUTPUT]
for target in targets:
    target.mkdir(parents=True,exist_ok=True)
    for f in files:
        destination=target/f.relative_to(HERE)
        destination.parent.mkdir(parents=True,exist_ok=True)
        if a.research and destination.exists() and destination.read_bytes()!=f.read_bytes():
            raise ValueError(f"different existing research file, not overwritten: {destination}")
        shutil.copyfile(f,destination)
    for entry in entries:
        actual=hashlib.sha256((target/entry["path"]).read_bytes()).hexdigest()
        if actual!=entry["sha256"]:raise ValueError(f"copy hash mismatch: {entry['path']}")
    if (target/"manifest.json").read_bytes()!=(HERE/"manifest.json").read_bytes():
        raise ValueError("manifest mismatch")
    print(json.dumps({"target":str(target),"verified_files":len(files),"bytes":sum(z["bytes"] for z in entries)}))
