"""E2b: Schwanz der Taylor-Koeffizienten von g(λ,x) um reelles λ0 (hohe Genauigkeit). Ausgabe |c_k| ρ^k und lokale Quotienten.
Aufruf: python e2b_schwanz.py DPS M lam0 rho x"""
import sys, json, time, mpmath as mp
from regulaer import g_lam
from motor import HIER
dps=int(sys.argv[1]); M=int(sys.argv[2]); mp.mp.dps=dps
l0=mp.mpf(sys.argv[3]); rho=mp.mpf(sys.argv[4]); x=mp.mpf(sys.argv[5])
t=time.time()
v=[g_lam(l0+rho*mp.expjpi(2*mp.mpf(j)/M),x) for j in range(M)]
c=[mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/M for k in range(M//2)]
print(f"lam0={l0} rho={rho} d=1-lam0={1-l0} t={time.time()-t:.0f}s")
prev=None
for k in range(0,M//2,8):
    a=abs(c[k])
    r = (abs(c[k])/abs(c[k-8]))**(mp.mpf(1)/8) if k>=8 else None
    print(k, mp.nstr(a,6), "" if r is None else f" rho/R≈{mp.nstr(r,5)}  R≈{mp.nstr(rho/r,5)}")
(HIER/"daten").mkdir(exist_ok=True)
(HIER/"daten"/f"e2b_{sys.argv[3]}_{sys.argv[4]}_{dps}_{M}.json").write_text(json.dumps([mp.nstr(abs(z),12) for z in c]))
