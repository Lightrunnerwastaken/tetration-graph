"""E3: Kneser knapp ausserhalb des Shell-Thron-Rands am goldenen Randpunkt.

t = (1+eps) e^{2 pi i alpha}, b = exp(t e^{-t}), oberer Fixpunkt L = e^t (Multiplikator t).
Schröder-Konstante kappa: a(T(z)-L) = v_t(kappa t^z) + O(e^{-2 pi Im z}); v_t = Poincaré-Funktion zu t.
Vorhersage (Siegel-Grenze): kappa -> lambda0 * phi(-lambda0) = lambda0^2 w*,  |kappa| -> r = 1.68854.
Aufruf: python e3_kneser_aussen.py eps [dps]
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
L = mp.exp(t)
bs = gpzahl(b, dps + 10)
zs = [mp.mpc(0.5, 3), mp.mpc(0.5, 4), mp.mpc(0.25, 3), mp.mpc(0.5, 10), mp.mpc(0.5, -3), mp.mpc(0.5, -10), mp.mpc(0.5, 0), mp.mpc(0, 0), mp.mpc(1, 0)]
ex = [f"sexp({gpzahl(z, dps + 5)})" for z in zs]
T = werte(bs, ex, dps)
schliessen()

# Poincaré-Funktion v_t
NC = 40
c = [mp.mpc(0)] * (NC + 1); e = [mp.mpc(0)] * (NC + 1)
c[1] = e[0] = e[1] = mp.mpc(1)
for n in range(2, NC + 1):
    R = mp.fsum(k * c[k] * e[n - k] for k in range(1, n)) / n
    c[n] = R / (t ** (n - 1) - 1)
    e[n] = c[n] + R
v = lambda w: mp.polyval(c[::-1], w)
dv = lambda w: mp.polyval([k * c[k] for k in range(NC, 0, -1)], w)

print(f"eps={eps}  b={mp.nstr(b, 15)}  L={mp.nstr(L, 15)}")
for z, Tz in zip(zs, T):
    print(f"  T({mp.nstr(z, 4)}) = {mp.nstr(Tz, 15)}   |T-L| = {mp.nstr(abs(Tz - L), 5)}")
for z, Tz in list(zip(zs, T))[:3]:
    target = a * (Tz - L)
    w = target
    for _ in range(50):
        w = w - (v(w) - target) / dv(w)
    kappa = w / t ** z
    print(f"  kappa aus z={mp.nstr(z, 4)}: {mp.nstr(kappa, 12)}   |kappa| = {mp.nstr(abs(kappa), 12)}")
