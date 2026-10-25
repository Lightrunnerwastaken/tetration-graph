"""E3: Ort der dominanten Singularität von θ_U aus den Koeffizienten (Rest r_n = a_n − b_n nach Abzug des q=1-Modells).
Quotient r_{n+1}/r_n -> ρ = 1/q*,  z* = log(q*)/(2πi) (mod 1);  Prony mit p Termen auf n in [n0, n1).
Aufruf: python e3_prony.py <welle-json> [n0 n1]
"""
import json
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 30
d = json.load(open(sys.argv[1], encoding="utf-8"))
n0 = int(sys.argv[2]) if len(sys.argv) > 2 else 120
n1 = int(sys.argv[3]) if len(sys.argv) > 3 else 200
w = d["wellen"]["U"]
a = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in w["a"]]
b = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in w["b"]]
r = [a[n] - b[n] for n in range(len(a))]
print(d["basis"], "a_0 =", mp.nstr(a[0], 10))
for n in [40, 80, 120, 160, n1 - 1]:
    q = r[n + 1] / r[n]
    z = mp.log(1 / q) / (2j * mp.pi)
    print(f"  n={n:4d}  r_n+1/r_n = {mp.nstr(q, 8)}  |.| = {mp.nstr(abs(q), 6)}  z* = {mp.nstr(mp.re(z) % 1, 5)} +({mp.nstr(mp.im(z), 5)})i   |r_n| = {mp.nstr(abs(r[n]), 4)}")
for p in [1, 2, 3]:
    M = n1 - n0
    X = np.array([[complex(r[n0 + i + j]) for j in range(p)] for i in range(M - p)])
    yv = -np.array([complex(r[n0 + i + p]) for i in range(M - p)])
    c, res, rk, sv = np.linalg.lstsq(X, yv, rcond=None)
    roots = np.roots(np.concatenate(([1], c[::-1])))
    zs = [np.log(1 / rt) / (2j * np.pi) for rt in roots]
    print(f"  Prony p={p}: z* = " + ", ".join(f"{z.real % 1:.5f}{z.imag:+.5f}i" for z in zs))
