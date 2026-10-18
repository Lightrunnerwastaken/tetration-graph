"""E1: Kneser T_b gegen die untere η-Tetration für b = η + δ ↓ η.  D(x,δ) = T_b(x) − T_η(x).
Aufruf: python e1_eta_grenze.py DPS d1 d2 ...   (d als GP-Ausdruck, z.B. 1/1000)"""
import sys, time, json, mpmath as mp
from motor import werte, schliessen, HIER
dps=int(sys.argv[1]); mp.mp.dps=dps
X=["0.5","-0.5","-1.5","2"]
out={}
for d in sys.argv[2:]:
    b=f"(exp(exp(-1))+{d})"
    t=time.time()
    v=werte(b,[f"sexp({x})" for x in X]+[f"sexpeta({x})" for x in X],dps)
    dd=mp.mpf(eval(d.replace('/','*mp.mpf(1)/'))) if '/' in d else mp.mpf(d)
    D=[mp.re(v[i]-v[i+len(X)]) for i in range(len(X))]
    print(d, " ".join(mp.nstr(Di/dd,14) for Di in D), f"{time.time()-t:.1f}s", flush=True)
    out[d]=[mp.nstr(Di,dps) for Di in D]
    schliessen()
(HIER/"daten").mkdir(exist_ok=True)
(HIER/"daten"/f"e1_{dps}.json").write_text(json.dumps({"X":X,"D":out},indent=1))
