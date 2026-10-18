"""E2: Konvergenzradius von λ ↦ g(λ,x) = T^reg(x) (Basis b = exp(λe^{−λ})) um Zentren λ0 im Einheitskreis.
Ist der Radius = Abstand zum Einheitskreis, so liegt dort eine Singularität. Aufruf: python e2_wurzeltest.py DPS M  re im rho_anteil x"""
import sys, mpmath as mp
from regulaer import g_lam
dps=int(sys.argv[1]); M=int(sys.argv[2]); mp.mp.dps=dps
l0=mp.mpc(mp.mpf(sys.argv[3]),mp.mpf(sys.argv[4])); q=mp.mpf(sys.argv[5]); x=mp.mpf(sys.argv[6])
d=1-abs(l0); rho=q*d
v=[g_lam(l0+rho*mp.expjpi(2*mp.mpf(j)/M),x) for j in range(M)]
c=[mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/M for k in range(M//2)]
print(f"λ0={mp.nstr(l0,6)} Abstand zum Kreis d={mp.nstr(d,6)}  rho={mp.nstr(rho,6)}")
for k in [8,16,24,32,40,48,56,64,80,96,112]:
    if k<M//2 and abs(c[k])>0:
        print(k, mp.nstr(abs(c[k]),5), " R_k/d =", mp.nstr(rho*abs(c[k])**(-mp.mpf(1)/k)/d,6))
