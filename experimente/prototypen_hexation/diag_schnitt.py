"""Diagnose: kreuzt beim Sprung von H(iy) (grobe θ-Lösung, y ≈ 1.4 → 1.6) ein Tetrations-Argument der Kette den Schnitt (−∞, −2]?"""
import sys, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from gemeinsam import kneser
from hexq import HexQ
from stufen import psi_wert
mp.mp.dps = 40
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
alt = json.load(open(HIER + "/prototypen_hexation/hex_lauf_2.json"))
g = [mp.mpf(v) for v in alt["32"]["g"]]
K = 32
def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]
gh = dft(g)
def g_at(x):
    return mp.re(gh[0] + sum(2 * mp.re(gh[n] * mp.expjpi(2 * n * x)) for n in range(1, K // 2)) + mp.re(gh[K // 2] * mp.cos(mp.pi * K * x)))
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
L, lam, c, psiT = HQ.L, HQ.lam, HQ.c, HQ.psiT
M = 8


def P_spur(ws):
    """Pentation wie stufen.pentation, aber mit Protokoll aller Tetrations-Argumente."""
    ys = [psi_wert(L, psiT, c * lam ** (w - M)) for w in ws]
    spur = [[] for _ in ws]
    for _ in range(M):
        for i, y in enumerate(ys):
            spur[i].append(y)
        ys = kneser("sexp", "2", ys)
    return ys, spur


def F_spur(zeta):
    s = mp.exp(zeta * HQ.lnmu)
    n = 0
    while abs(s / HQ.mu ** n) > HQ.rlok:
        n += 1
    y = psi_wert(HQ.q, HQ.psi, s / HQ.mu ** n)
    alle = []
    for k in range(n):
        (y,), (sp,) = P_spur([y])
        alle.append((k, sp))
    return y, alle


for yy in [mp.mpf("1.30") + mp.mpf(k) / 50 for k in range(0, 21)]:
    zeta = XS + 1j * yy + theta(XS + 1j * yy)
    w, alle = F_spur(zeta)
    nahe = [(k, j, mp.nstr(t, 5)) for k, sp in alle for j, t in enumerate(sp) if mp.re(t) < -1.9 and abs(mp.im(t)) < 0.6]
    print(f"y = {mp.nstr(yy, 3)}  H = {mp.nstr(w, 8)}  Argumente nahe am Schnitt: {nahe}", flush=True)
