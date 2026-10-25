"""Gebiet 26: Raabe-Binet-Form auf jeder Geraden und erstes senkrechtes Moment für Basis 2 (Gegenprobe der allgemeinen Formel).

Allgemeine Basis b > η, λ = L ln b, D(t) = b^s − s, s = sexp(x0+it), K(x) = ∫_0^∞ (sexp(x+it) − L) dt, Λ(x) = ∫_x^{x+1} ln sexp':
  RBb  i ln b · K(x) + Λ(x) − x log λ  ≡  ln κ_b + ½ log λ + log log λ   (mod 2πi), unabhängig von x
  M1b  ln b · Im ∫_0^∞ t D dt − ln b · Im L (x0 + 1)  ≡  arg κ_b + arg log λ   (mod 2π)
  Vb   ln b · Re K(x0) − ln b · Im L (x0 + ½)        ≡  arg κ_b + arg log λ   (mod 2π)
Daten nur lesend: etappe3/daten/taylor_b2_80_256.json, etappe2/daten/a_2_80_0.003_16384.json.
Aufruf: python s2_basis2.py [DPS]
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E3 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe3"
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = dps
b = mp.mpf(2)
lb = mp.log(b)
L = -mp.lambertw(-lb) / lb
if mp.im(L) < 0:
    L = mp.conj(L)
lam = L * lb
llam = mp.log(lam)
K = 70
KAPPA2_E3 = mp.mpc("-0.0912950128173059192138757198914", "-3.00920981975562728783175111341")

ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in json.loads((E3 / "daten" / "taylor_b2_80_256.json").read_text(encoding="utf-8")).items()}
DZEN = {c: [k * a for k, a in enumerate(v)][1:] for c, v in ZEN.items()}

P = [L, mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
for k in range(2, K + 1):
    q = [mp.mpc(0)] + [lb * x for x in P[1:k]] + [mp.mpc(0)] * (K + 1 - k)
    e = [mp.mpc(1)] + [mp.mpc(0)] * k
    for m in range(1, k + 1):
        e[m] = mp.fsum(j * q[j] * e[m - j] for j in range(1, m + 1)) / m
    P[k] = L * e[k] / (lam ** k - lam)


def Psi(u):
    n = 0
    while abs(u) > mp.mpf("0.02"):
        u /= lam
        n += 1
    w = mp.polyval(P[::-1], u)
    for _ in range(n):
        w = b ** w
    return w


d = json.loads((E3.parent / "etappe2" / "daten" / "a_2_80_0.003_16384.json").read_text(encoding="utf-8"))
A = [mp.mpc(mp.mpf(x[0]), mp.mpf(x[1])) for x in d["a"]]
_memo = {}


def sexp(x0, t):
    key = (mp.nstr(x0, dps + 5), mp.nstr(t, dps + 5))
    if key in _memo:
        return _memo[key]
    z = mp.mpc(x0, t)
    w = None
    if t <= mp.mpf("0.09"):
        c = min(ZEN, key=lambda c: abs(z - c))
        if abs(z - c) <= mp.mpf("0.1"):
            w = mp.polyval(ZEN[c][::-1], z - c)
    if w is None:
        qq = mp.expjpi(2 * z); th = mp.mpc(0); qn = mp.mpc(1)
        for ak in A:
            th += ak * qn; qn *= qq
            if abs(qn) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        w = Psi(mp.exp(llam * (z + th)))
    _memo[key] = w
    return w


def reell(x):
    x = mp.mpf(x)
    if x <= mp.mpf("0.5"):
        c = min(ZEN, key=lambda c: abs(x - c))
        return mp.polyval(ZEN[c][::-1], x - c), mp.log(mp.polyval(DZEN[c][::-1], x - c))
    v1, l1 = reell(x - 1)
    return b ** v1, mp.log(lb) + v1 * lb + l1


KN = [0, mp.mpf("0.09"), mp.mpf("0.5"), 2, 6, 15, 35, 70, 140]


def mod2pi(x):
    return x - 2 * mp.pi * mp.nint(x / (2 * mp.pi))


if __name__ == "__main__":
    print(f"Basis 2, dps={dps}, L = {mp.nstr(L, 25)}, λ = {mp.nstr(lam, 25)}")
    soll_c = mp.log(KAPPA2_E3) + llam / 2 + mp.log(llam)
    soll_a = mp.arg(KAPPA2_E3) + mp.im(mp.log(llam))
    consts = []
    for x0 in (mp.mpf(0), mp.mpf("-0.4"), mp.mpf("0.2")):
        Kx = mp.quad(lambda t: sexp(x0, t) - L, KN)
        W1 = mp.quad(lambda t: t * (b ** sexp(x0, t) - sexp(x0, t)), KN)
        pts = [x0] + [p for p in (mp.mpf("-0.3"), mp.mpf("-0.1"), mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.5"), mp.mpf("0.7"), mp.mpf("0.9"), mp.mpf("1.1")) if x0 < p < x0 + 1] + [x0 + 1]
        Lam = mp.quad(lambda x: reell(x)[1], pts)
        c = 1j * lb * Kx + Lam - x0 * llam
        consts.append(c)
        dc = c - soll_c
        k = mp.nint(mp.im(dc) / (2 * mp.pi))
        print(f"x0={mp.nstr(x0,2):5s} RBb: Konstante = {mp.nstr(c, 28)}; minus Soll = {mp.nstr(dc, 5)} → Zweigzahl {k}, Rest {mp.nstr(abs(dc - 2j * mp.pi * k), 3)}")
        m1 = lb * mp.im(W1) - lb * mp.im(L) * (x0 + 1)
        v = lb * mp.re(Kx) - lb * mp.im(L) * (x0 + mp.mpf(1) / 2)
        print(f"          M1b: {mp.nstr(m1, 28)} − (arg κ₂ + arg log λ) = {mp.nstr(m1 - soll_a, 5)} (mod 2π: {mp.nstr(mod2pi(m1 - soll_a), 3)})")
        print(f"          Vb : {mp.nstr(v, 28)} − (arg κ₂ + arg log λ) = {mp.nstr(v - soll_a, 5)} (mod 2π: {mp.nstr(mod2pi(v - soll_a), 3)})", flush=True)
    print(f"Streuung der RBb-Konstante über drei Geraden: {mp.nstr(max(abs(a - c) for a in consts for c in consts), 3)}")
