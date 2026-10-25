"""E5: Prüft, ob die Motorwerte die Gleichung T(z+1) = b^T(z) erfüllen (abseits ganzer Höhen).
Aufruf: python e5_fgl_probe.py eps [dps]
"""
import sys
import mpmath as mp
from motor import werte, schliessen, gpzahl

eps = mp.mpf(sys.argv[1])
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = dps + 15
alpha = (3 - mp.sqrt(5)) / 2
t = (1 + eps) * mp.exp(2j * mp.pi * alpha)
a = t * mp.exp(-t)
b = mp.exp(a)
zs = [mp.mpc(0.5, 0), mp.mpc(0.3, 0.5), mp.mpc(0.5, -1), mp.mpc(-0.5, 0.2)]
ex = [f"sexp({gpzahl(z, dps + 5)})" for z in zs] + [f"sexp({gpzahl(z + 1, dps + 5)})" for z in zs]
T = werte(gpzahl(b, dps + 10), ex, dps)
schliessen()
n = len(zs)
for i, z in enumerate(zs):
    print(f"eps={mp.nstr(eps,4)} z={mp.nstr(z,3)}: |T(z+1) - b^T(z)| = {mp.nstr(abs(T[n+i] - mp.exp(a*T[i])), 3)}")
