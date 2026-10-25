"""Rang 6: Kneser-Hexation (Motor, über eta5) auf [-1,0] als Tschebyschow-Interpolation; Basiswechsel-Phase Phi6.

Hex(T)                hex_b auf der reellen Achse: hex(x+1) = P_b(hex(x)), hex(0) = 1; ihex = Umkehrung.
G5x(Tb, Td, w)        = spog_b(P_d(w)) exakt (auch für sehr grosse w).
H6(Tb, Td, y)         = lim spog_b^n(P_d^n(y));  Phi6(th) = ihex_b(H6(hex_d(th))) - th.
"""
import json

import mpmath as mp

import lib


class Hex:
    def __init__(self, T, N=120, dps=60):
        d = json.loads((lib.CDIR / f"hexx_{lib._tag(T.name)}_{N}_{dps}.json").read_text(encoding="utf-8"))
        self.T = T
        f = [mp.mpf(v) for v in d["T"]]
        c = []
        for k in range(N):
            s = mp.fsum(f[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / N) for j in range(N))
            c.append(2 * s / N)
        c[0] /= 2
        self.c = c
        self.dc = lib.Tet._cheb_abl(c)

    def h0(self, x):
        return lib.Tet._clen(self.c, 2 * x + 1)

    def dh0(self, x):
        return 2 * lib.Tet._clen(self.dc, 2 * x + 1)

    def hex(self, x):
        m = int(mp.ceil(x))
        y = self.h0(x - m)
        for _ in range(m):
            y = self.T.P(y)
        for _ in range(-m):
            y = self.T.spog(y)
        return y

    def ihex(self, y):
        n = 0
        while y > 1:
            y = self.T.spog(y)
            n += 1
        while y <= 0:
            y = self.T.P(y)
            n -= 1
        t = y - 1
        for _ in range(200):
            dt = (self.h0(t) - y) / self.dh0(t)
            t -= dt
            if abs(dt) < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        return n + t


def Hn(Tb, Td, y, n, X=None):
    """slog_b^n(T_d^n(y)) exakt; für grosse n gleich dem Limes H(y)."""
    if X is None:
        X = mp.mp.dps // 2
    if n == 0:
        return y
    A = Td.slog(mp.mpf(10) ** X)
    U = [y]
    while len(U) < n and U[-1] < A:
        U.append(Td.T(U[-1]))
    w = lib.R(Tb, Td, U[-1])
    for _ in range(len(U) - 1):
        w = Tb.slog(w)
    return w


def G5x(Tb, Td, w):
    n = int(mp.floor(w))
    th = w - n
    if n < 0:
        return Tb.spog(Td.P(w))
    return n + Tb.spog(Hn(Tb, Td, Td.P(th), n))


def H6(Tb, Td, y, X=None):
    if X is None:
        X = mp.mp.dps // 3
    A = Td.spog(mp.mpf(10) ** X)
    V = [y]
    while V[-1] < A:
        V.append(Td.P(V[-1]))
        if len(V) > 40:
            raise RuntimeError("H6: Kette zu lang")
    w = G5x(Tb, Td, V[-1])
    for _ in range(len(V) - 1):
        w = Tb.spog(w)
    return w


def phi6(Tb, Td, Hb, Hd, th, X=None):
    return Hb.ihex(H6(Tb, Td, Hd.hex(th), X)) - th


def G6_endlich(Tb, Td, Hb, Hd, th, n):
    return Hb.ihex(Hd.hex(n + th)) - n - th
