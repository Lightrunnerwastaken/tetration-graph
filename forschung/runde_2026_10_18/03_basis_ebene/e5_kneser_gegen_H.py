"""E5: Kneser T_{η+δ}(x) (Motor) gegen die Taylor-Reihe H_N(δ)=Σ a_n δ^n bei η aus der regulären Familie (e4-Datei).
Aufruf: python e5_kneser_gegen_H.py e4datei x NGUT d1 d2 ..."""
import sys, json, mpmath as mp
from motor import werte, schliessen, HIER
mp.mp.dps=50
a=[mp.mpf(z) for z in json.loads((HIER/"daten"/sys.argv[1]).read_text())]
x=sys.argv[2]; N=int(sys.argv[3]); import os; KDPS=int(os.environ.get("KDPS","40"))
for d in sys.argv[4:]:
    D=mp.mpf(d.split('/')[0])/mp.mpf(d.split('/')[1])
    K=mp.re(werte(f"(exp(exp(-1))+{d})",[f"sexp({x})"],int(sys.argv[-1]) if False else KDPS)[0])
    H=mp.fsum(a[n]*D**n for n in range(N))
    print(f"δ={d}: Kneser={mp.nstr(K,22)}  K−H_{N}={mp.nstr(K-H,4)}  letztes Glied={mp.nstr(a[N-1]*D**(N-1),3)}", flush=True)
    schliessen()
