"""Small selected critical-family samples; each GP call has a timeout.

No new critical roots, no interval claims, no shared caches. Reuses the frozen
preceding-round engine and stores runtime inputs and raw values in this round.
"""
import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).absolute().parent
RUNTIME = HERE / "runtime"
sys.path.insert(0,str(RUNTIME))
from common import engine, mp, stringify

p=argparse.ArgumentParser()
p.add_argument("--rank",type=int,required=True)
p.add_argument("--dps",type=int,default=20)
p.add_argument("--max-seconds",type=int,default=240)
a=p.parse_args()
mp.mp.dps = a.dps+30
with (HERE/"inputs/eta_values.csv").open(encoding="utf-8-sig",newline="") as f:
    row=next(z for z in csv.DictReader(f) if int(z["map_rank"])==a.rank)
start=time.monotonic()
gp=engine(a.dps,RUNTIME/"hyper.gp",logged=True,timeout=a.max_seconds-30)
out={"rank":a.rank,"requested_dps":a.dps,"base":row["eta"],"fold_height":row["e"],
     "base_source_sha256":hashlib.sha256((HERE/"inputs/eta_values.csv").read_bytes()).hexdigest(),
     "status":"running","values":[]}
target=HERE/f"critical_profile_r{a.rank}_{a.dps}.json"
def save():
    out["seconds"]=time.monotonic()-start
    target.write_text(json.dumps(out,indent=2),encoding="utf-8")
try:
    gp.eval_batch(row["eta"],["hrdebug=1", f"hrmaxdigits[7]={a.dps-6}",f"hinit({a.rank})"])
    out["engine_metadata"]={"expressions":["precis",f"hrerror[{a.rank}]",f"hrkind[{a.rank}]"]}
    out["engine_metadata"]["values"]=[stringify(z,a.dps+10) for z in gp.eval_batch(row["eta"],out["engine_metadata"]["expressions"])]
    print("initialized",a.rank,time.monotonic()-start,flush=True)
    xs=["-.5","-.25","0",".25",".5","1","1.5","2","3",row["e"],
        *[str(-k) for k in range(1,a.rank-2)]]
    for x in xs:
        if time.monotonic()-start > a.max_seconds-20:
            raise TimeoutError("profile budget reached")
        expr=f"hval({a.rank},{x})"
        y=gp.eval_batch(row["eta"],[expr])[0]
        out["values"].append({"x":x,"expression":expr,"value":stringify(y,a.dps+10)})
        save()
        print(a.rank,x,stringify(y,17),flush=True)
    # An independent recursion evaluation tests samples, including left side.
    exprs=[f"hval({a.rank},{x}+1)-hval({a.rank-1},hval({a.rank},{x}))" for x in ["-.5",".5"]]
    out["recursion_checks"]={"expressions":exprs,"values":[stringify(z,a.dps+10) for z in gp.eval_batch(row["eta"],exprs)]}
    out["status"]="completed"
except Exception as exc:
    out["status"]="failed_or_partial";out["error"]=f"{type(exc).__name__}: {exc}"
    print(out["error"],flush=True)
finally:
    gp.close();save()
    print("seconds",out["seconds"],flush=True)
