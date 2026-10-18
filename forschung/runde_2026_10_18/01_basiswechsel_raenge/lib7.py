"""Rang 7: hyper(7, x) aus dem Motor (regulär am linken reellen Fixpunkt der Hexation) auf [-1,0] als Tschebyschow-Interpolation;
Basiswechsel-Phase Phi7(th) = lim A7_b(H7_d(n+th)) - n - th.
"""
import json

import mpmath as mp

import lib
import lib6


class Hyp7:
    def __init__(self, Hx, N=120, dps=60):
        d = json.loads((lib.CDIR / f"hyper7_{lib._tag(Hx.T.name)}_{N}_{dps}.json").read_text(encoding="utf-8"))
        self.Hx = Hx
        f = [mp.mpf(v) for v in d["T"]]
        c = []
        for k in range(N):
            s = mp.fsum(f[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / N) for j in range(N))
            c.append(2 * s / N)
        c[0] /= 2
        self.c = c
        self.dc = lib.Tet._cheb_abl(c)
        self.schwanz = max(abs(x) for x in c[-10:])

    def h0(self, x):
        return lib.Tet._clen(self.c, 2 * x + 1)

    def dh0(self, x):
        return 2 * lib.Tet._clen(self.dc, 2 * x + 1)

    def hyp(self, x):
        m = int(mp.ceil(x))
        y = self.h0(x - m)
        for _ in range(m):
            y = self.Hx.hex(y)
        for _ in range(-m):
            y = self.Hx.ihex(y)
        return y

    def ihyp(self, y):
        n = 0
        while y > 1:
            y = self.Hx.ihex(y)
            n += 1
        while y <= 0:
            y = self.Hx.hex(y)
            n -= 1
        t = y - 1
        for _ in range(200):
            dt = (self.h0(t) - y) / self.dh0(t)
            t -= dt
            if abs(dt) < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        return n + t


def H6n(Tb, Td, y, n, X=None):
    """spog_b^n(P_d^n(y)) exakt; für grosse n gleich dem Limes H6(y)."""
    if X is None:
        X = mp.mp.dps // 3
    if n == 0:
        return y
    A = Td.spog(mp.mpf(10) ** X)
    V = [y]
    while len(V) < n and V[-1] < A:
        V.append(Td.P(V[-1]))
    w = lib6.G5x(Tb, Td, V[-1])
    for _ in range(len(V) - 1):
        w = Tb.spog(w)
    return w


def G6x(Tb, Td, Hb, Hd, w):
    """ihex_b(hex_d(w)) exakt, auch für sehr grosse w."""
    n = int(mp.floor(w))
    th = w - n
    if n < 0:
        return Hb.ihex(Hd.hex(w))
    return n + Hb.ihex(H6n(Tb, Td, Hd.hex(th), n))


def H7(Tb, Td, Hb, Hd, y, X=None):
    if X is None:
        X = mp.mp.dps // 4
    A = Hd.ihex(mp.mpf(10) ** X)
    Z = [y]
    while Z[-1] < A:
        Z.append(Hd.hex(Z[-1]))
        if len(Z) > 40:
            raise RuntimeError("H7: Kette zu lang")
    w = G6x(Tb, Td, Hb, Hd, Z[-1])
    for _ in range(len(Z) - 1):
        w = Hb.ihex(w)
    return w


def phi7(Tb, Td, Hb, Hd, Yb, Yd, th, X=None):
    return Yb.ihyp(H7(Tb, Td, Hb, Hd, Yd.hyp(th), X)) - th
