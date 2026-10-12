"""Freeze user-supplied research material and previous numerical data."""
import hashlib
import json
import shutil
from pathlib import Path

here=Path(__file__).absolute().parent
inputs=here/"inputs"
inputs.mkdir(exist_ok=True)
sources=[
    Path(r"D:\Users\janis\Downloads\hyperoperation_eta_2_4_successor.md"),
    Path(r"D:\Users\janis\Downloads\hyperoperation_high_rank_convergence_eta_e.md"),
    here.parents[1]/"outputs/hyper_eta/eta_values.csv",
]
records=[]
for source in sources:
    target=inputs/source.name
    shutil.copyfile(source,target)
    records.append({"original":str(source),"frozen_copy":"inputs/"+source.name,
                    "sha256":hashlib.sha256(target.read_bytes()).hexdigest()})
(here/"input_provenance.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
print(json.dumps(records,indent=2))
