"""E7: Schwanz der b-Taylor-Koeffizienten der regulären Familie um b0 = η − s: P_k = |c_k| s^k (Radius s ⇒ P_k konstant).
Aufruf: python e7_schwanz_b.py DPS M q x s1 s2 ..."""
import sys, json, mpmath as mp
from regulaer import T_reg
from motor import HIER
dps=int(sys.argv[1]); M=int(sys.argv[2]); mp.mp.dps=dps
q=mp.mpf(sys.argv[3]); x=mp.mpf(sys.argv[4]); eta=mp.exp(mp.exp(-1)); res={}
for sv in sys.argv[5:]:
    s=mp.mpf(sv); b0=eta-s; rho=q*s
    v=[T_reg(b0+rho*mp.expjpi(2*mp.mpf(j)/M),x) for j in range(M)]
    c=[mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/M*(s/rho)**k for k in range(M//2)]
    P=[abs(z) for z in c]
    ks=list(range(8,M//2,8))
    print(f"s={sv}:", " ".join(f"{k}:{mp.nstr(P[k],3)}" for k in ks), flush=True)
    res[sv]=[mp.nstr(z,8) for z in P]
(HIER/"daten").mkdir(exist_ok=True)
(HIER/"daten"/f"e7_{sys.argv[3]}_{sys.argv[4]}_{dps}_{M}.json").write_text(json.dumps(res))
