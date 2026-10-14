"""Sonde: θ(z) = h(sexp z) − z nahe der reellen Achse, Basis e.
(1) senkrecht auf z = iy, y = 10^-k  (2) waagrecht z = x + i·y0.
Aufruf: python p1_sonde.py <dps>
"""
import sys
import mpmath as mp
from kern import Fix, sexp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 50
mp.mp.dps = DPS + 20
fx = Fix(mp.e, K=40)
L = fx.L
per = 2j * mp.pi / mp.log(fx.lam)
print("L =", mp.nstr(L, 20), " Periode 2πi/L =", mp.nstr(per, 12))

ys = [mp.mpf(10) ** (-mp.mpf(k) / 4) for k in range(2, 4 * 30 + 1)]
zs = [1j * y for y in ys]
ws = sexp("exp(1)", zs, DPS)
prev = None
for y, z, w in zip(ys, zs, ws):
    th = fx.h(w) - z
    if prev is not None:
        m = min(range(-5, 6), key=lambda m: abs(th + m * per - prev))
        th += m * per
    prev = th
    print(f"y=1e{float(mp.log10(y)):7.2f}  w-1={mp.nstr(w-1, 6):>28}  theta={mp.nstr(th, 12)}")
