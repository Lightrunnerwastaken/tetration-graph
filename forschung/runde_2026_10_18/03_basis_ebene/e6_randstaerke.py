"""E6: Stärke der Randsingularität von g(λ,x) entlang |λ|=1: Zentren λ0 = r0 e^{iθ}, Kreis ρ; geschätzter Schwanz-Radius R und
Amplitude P = |c_k| ρ^{-k}... genauer P_k = |c_k| d^k (c_k Koeff. in (λ−λ0)), d = 1−r0.  Aufruf: python e6_randstaerke.py DPS M r0 rho x th1 th2 ...  (θ in Einheiten π)"""
import sys, json, mpmath as mp
from regulaer import g_lam
from motor import HIER
dps=int(sys.argv[1]); M=int(sys.argv[2]); mp.mp.dps=dps
r0=mp.mpf(sys.argv[3]); rho=mp.mpf(sys.argv[4]); x=mp.mpf(sys.argv[5])
d=1-r0; res={}
for th in sys.argv[6:]:
    l0=r0*mp.expjpi(mp.mpf(th))
    v=[g_lam(l0+rho*mp.expjpi(2*mp.mpf(j)/M),x) for j in range(M)]
    c=[mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/M for k in range(M//2)]  # = koeff * rho^k
    P=[abs(c[k])*(d/rho)**k for k in range(M//2)]
    ks=list(range(M//4,M//2,8))
    R=[rho/((abs(c[k])/abs(c[k-16]))**(mp.mpf(1)/16)) for k in ks]
    print(f"θ={th}π: R/d an k={ks}:", " ".join(mp.nstr(r/d,3) for r in R), "| P_k:", " ".join(mp.nstr(P[k],3) for k in ks), flush=True)
    res[th]=[mp.nstr(z,8) for z in P]
(HIER/"daten").mkdir(exist_ok=True)
(HIER/"daten"/f"e6_{sys.argv[3]}_{sys.argv[4]}_{dps}_{M}.json").write_text(json.dumps(res))
