"""E6: Schröder-Konstante kappa(eps) der Motor-Lösung (Kneser-Art) entlang t = (1+eps) e^{2 pi i alpha}, goldener Randpunkt.
Zusätzlich reguläre Konstante (Koenigs) für eps<0.  Ausgabe: cache/e6_scan.txt
Aufruf: python e6_scan.py eps1 eps2 ...
"""
import sys
import mpmath as mp
from motor import werte, schliessen, gpzahl

dps = 30
mp.mp.dps = 45
alpha = (3 - mp.sqrt(5)) / 2


def poincare(t, NC=40):
    c = [mp.mpc(0)] * (NC + 1); e = [mp.mpc(0)] * (NC + 1)
    c[1] = e[0] = e[1] = mp.mpc(1)
    for n in range(2, NC + 1):
        R = mp.fsum(k * c[k] * e[n - k] for k in range(1, n)) / n
        c[n] = R / (t ** (n - 1) - 1)
        e[n] = c[n] + R
    v = lambda w: mp.polyval(c[::-1], w)
    dv = lambda w: mp.polyval([k * c[k] for k in range(NC, 0, -1)], w)
    return v, dv


out = open("cache/e6_scan.txt", "a", encoding="utf-8")
for s in sys.argv[1:]:
    eps = mp.mpf(s)
    t = (1 + eps) * mp.exp(2j * mp.pi * alpha)
    a = t * mp.exp(-t)
    b = mp.exp(a)
    L = mp.exp(t)
    z = mp.mpc(0.5, 3)
    T3, T0 = werte(gpzahl(b, dps + 10), [f"sexp({gpzahl(z, 35)})", "sexp(0.5)"], dps)
    schliessen()
    v, dv = poincare(t)
    target = a * (T3 - L)
    w = target
    for _ in range(40):
        w = w - (v(w) - target) / dv(w)
    kap = w / t ** z
    line = f"{float(eps):+.4f}  kappaK = {mp.nstr(kap, 10)}  |kappaK| = {mp.nstr(abs(kap), 10)}  T(0.5) = {mp.nstr(T0, 10)}"
    print(line, flush=True)
    out.write(line + "\n")
out.close()
