"""Kleine Potenzreihen-Werkzeuge für Strang P (mpmath). Reihen = Listen von Koeffizienten (Länge K+1)."""
import json
from pathlib import Path

import mpmath as mp

GEM = Path(__file__).resolve().parent.parent / "gemeinsam"


def kneser_koeff(n=4000):
    return [mp.mpf(s) for s in json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text())["a"]][:n + 1]


def fixpunkt_L():
    return mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))


def mul(a, b, K):
    out = [mp.mpf(0)] * (K + 1)
    for i, x in enumerate(a[:K + 1]):
        if x == 0:
            continue
        for j in range(0, K + 1 - i):
            if j < len(b):
                out[i + j] += x * b[j]
    return out


def add(a, b, K, s=1):
    out = []
    for k in range(K + 1):
        x = a[k] if k < len(a) else 0
        y = b[k] if k < len(b) else 0
        out.append(x + s * y)
    return out


def scal(c, a):
    return [c * x for x in a]


def inv(a, K):
    b = [1 / a[0]]
    for k in range(1, K + 1):
        s = 0
        for j in range(1, min(k, len(a) - 1) + 1):
            s += a[j] * b[k - j]
        b.append(-s / a[0])
    return b


def div(a, b, K):
    return mul(a, inv(b, K), K)


def deriv(a):
    return [k * a[k] for k in range(1, len(a))] or [mp.mpf(0)]


def shift(c, z0, K):
    """Taylor-Koeffizienten von sum c_k x^k an der Stelle z0, Ordnung <= K (Horner-artig über Binomialsumme)."""
    n = len(c)
    out = []
    # b_j = sum_k binom(k,j) c_k z0^(k-j)
    pw = [mp.mpf(1)]
    for _ in range(n):
        pw.append(pw[-1] * z0)
    for j in range(K + 1):
        s = 0
        bk = mp.mpf(1)  # binom(k,j) for k=j
        for k in range(j, n):
            s += bk * c[k] * pw[k - j]
            bk = bk * (k + 1) / (k + 1 - j)
        out.append(s)
    return out


def compose(outer_at_c, inner, K):
    """outer(inner(x)) als Reihe bis x^K; outer_at_c = Taylor-Koeffizienten von outer an c = inner(0)."""
    y = [mp.mpf(0)] + list(inner[1:K + 1])
    res = [outer_at_c[-1]]
    for co in reversed(outer_at_c[:-1]):
        res = add(mul(res, y, K), [co], K)
    return (res + [mp.mpf(0)] * (K + 1))[:K + 1]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def solve(A, b):
    return mp.lu_solve(A, b)


def exp_reihe(K):
    return [1 / mp.factorial(k) for k in range(K + 1)]
