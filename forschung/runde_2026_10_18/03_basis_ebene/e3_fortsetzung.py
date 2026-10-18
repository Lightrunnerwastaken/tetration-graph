"""E3: Taylor-Reihe der regulären Tetration b ↦ T^reg_b(x) um b0 = η − s (Cauchy-Kreis Radius q·s),
über η hinaus ausgewertet und mit Kneser (Motor) bei b = η + δ verglichen.
Aufruf: python e3_fortsetzung.py DPS M s q x  d1 d2 ..."""
import sys, json, mpmath as mp
from regulaer import T_reg
from motor import werte, schliessen, HIER
dps=int(sys.argv[1]); M=int(sys.argv[2]); mp.mp.dps=dps
s=mp.mpf(sys.argv[3]); q=mp.mpf(sys.argv[4]); x=sys.argv[5]; ds=sys.argv[6:]
eta=mp.exp(mp.exp(-1)); b0=eta-s; rho=q*s
v=[T_reg(b0+rho*mp.expjpi(2*mp.mpf(j)/M),mp.mpf(x)) for j in range(M)]
c=[mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/M/rho**k for k in range(M//2)]
print("Koeffizienten |c_k| s^k (Radius in Einheiten s):")
for k in [4,8,16,24,32,40,48,56,63]:
    if k<M//2: print(" ",k, mp.nstr(abs(c[k])*s**k,5), mp.nstr((abs(c[k]))**(-mp.mpf(1)/k)/s,5))
for d in ds:
    D=mp.mpf(eval(d)) if '/' not in d else mp.mpf(d.split('/')[0])/mp.mpf(d.split('/')[1])
    K=mp.re(werte(f"(exp(exp(-1))+{d})",[f"sexp({x})"],min(dps,40))[0])
    t=D+s
    part=[mp.fsum(c[k]*t**k for k in range(N)) for N in (16,24,32,40,48,56,M//2)]
    print("δ=",d," Kneser",mp.nstr(K,25))
    for N,p in zip((16,24,32,40,48,56,M//2),part):
        print("   N=",N," Reihe−Kneser =",mp.nstr(mp.re(p)-K,5), " Im",mp.nstr(mp.im(p),3))
schliessen()
