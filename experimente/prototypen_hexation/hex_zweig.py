"""Verlauf von H(iy) (grobe Lösung) und P^{-1}(H(iy)) entlang y — Sprünge des Hauptzweigs?"""
import sys, json
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/experimente")
import mpmath as mp
from hexq import HexQ
from hexation import Pentation
mp.mp.dps = 40
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
Pe = Pentation("2")
alt = json.load(open(r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad/hex_lauf_2.json"))
g = [mp.mpf(v) for v in alt["32"]["g"]]
K = 32
def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]
gh = dft(g)
def g_at(x):
    s = gh[0] + sum(2 * mp.re(gh[n] * mp.expjpi(2 * n * x)) for n in range(1, K // 2)) + mp.re(gh[K // 2] * mp.cos(mp.pi * K * x))
    return mp.re(s)
xs = [mp.mpf(j) / K for j in range(K)]
u = [mp.mpf(0)] * K
for it in range(200):
    v = [g_at(x + uu) for x, uu in zip(xs, u)]
    vh = dft(v)
    a = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
    un = [mp.re(mp.fsum(a[n] * mp.expjpi(2 * n * x) for n in range(len(a)))) for x in xs]
    if max(abs(p - r) for p, r in zip(u, un)) < mp.mpf(10) ** -22:
        break
    u = un
theta = lambda z: mp.fsum(a[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(a)))
XS = mp.mpf(alt["32"]["x_stern"])
ys = [mp.mpf(k) / 5 for k in range(0, 41)]
Hy = HQ.F([XS + 1j * yy + theta(XS + 1j * yy) for yy in ys])
Iy = Pe.inv(Hy)
Py = Pe.P(Hy)
for yy, h, i, p in zip(ys, Hy, Iy, Py):
    print(f"y = {mp.nstr(yy, 3):>4}  H = {mp.nstr(h, 8):>32}  P^-1(H) = {mp.nstr(i, 8):>32}  P(H) = {mp.nstr(p, 8)}")
