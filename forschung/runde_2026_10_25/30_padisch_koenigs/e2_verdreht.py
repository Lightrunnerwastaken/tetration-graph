"""E2: geschlossene Form fuer beliebige a (p teilt a nicht), auch a ≢ 1 mod p.
Muster zeta_n = omega(a)^(x_n mod (p-1)) / u  (x_n = a^^n, n >= -1), u = omega(a)^c.
   y_{-1} = -lam,  y_{n+1} = F_{zeta_n}(y_n) = lam (zeta_n e^{y_n} - 1),  n0 = ab hier zeta = 1,
   kappa = sigma_lam(y_{n0}) / (L lam^{n0}),   lam = u L Y(u L)  (Lambert),  alpha = u Y(u L).
Vergleich mit dem direkten Grenzwert (a^^n - alpha)/lam^n ueber exakte Tuerme mod p^N (wie Bericht 04, E2c).
"""
import sys, time
from math import factorial
from fractions import Fraction as Fr
from kr_lib import *
import ar_lib

t0 = time.time()
N = 60          # Arbeitsgenauigkeit fuer die Reihen
KM = 70         # Grad in y und lam
S = sk_series(KM, KM)
EXP = [Fr(1, factorial(k)) for k in range(KM + 1)]


class P:
    """p-adische Zahl als R/p^E mod p^N (absolut)."""
    pass


def to_int(R, E, p, Nw):
    # nur fuer E == 0 benutzen
    assert E == 0
    return R % p ** Nw


def padic_exp(y, p, Nw):
    R, E = eval_padic(EXP, y, p, Nw)
    return to_int(R, E, p, Nw)


def sigma_val(lam, y, p, Nw):
    """sigma_lam(y) mod p^Nw, v(y) >= 1; Rueckgabe ganze Zahl (Wert ist ganz, da |s_k| <= |1/k!|)."""
    q = p ** Nw
    tot = 0
    for k in range(1, KM + 1):
        Rk, Ek = eval_padic(S[k], lam, p, Nw + 30)   # s_k(lam) = Rk / p^Ek
        yk = pow(y, k, p ** (Nw + 60))
        t = Rk * yk
        assert t % p ** Ek == 0
        tot = (tot + t // p ** Ek) % q
    return tot


def kappa_closed(a, p, Nw):
    q = p ** Nw
    W = Nw + 20
    L = ar_lib.padic_log(a, p, W)
    if p == 2:
        u = 1 if a % 4 == 1 else -1
        om = u
        res = lambda x: x % 2
    else:
        om = ar_lib.teichmueller(a, p, W)
        c = ar_lib.tower_mod(a, 60, p - 1)
        u = pow(om, c, p ** W)
        res = lambda x: x % (p - 1)
    uinv = pow(u % p ** W, -1, p ** W)
    # Muster: x_n mod (p-1) fuer n = -1, 0, 1, ...
    xs = [0, 1] + [ar_lib.tower_mod(a, n, p - 1 if p > 2 else 2) for n in range(1, 40)]
    zetas = [pow(om % p ** W, res(x), p ** W) * uinv % p ** W for x in xs]  # zeta_{-1}, zeta_0, ...
    n0i = max([i for i, z in enumerate(zetas) if z != 1], default=-1) + 1   # Index in Liste
    n0 = n0i - 1                                                           # Hoehe
    Y = lambert_series(KM)
    uL = u * L % p ** W
    R, E = eval_padic(Y, uL, p, W)
    alpha = u * R % p ** W
    lam = L * alpha % p ** W
    y = (-lam) % p ** W
    for i in range(n0i):
        y = lam * (zetas[i] * padic_exp(y, p, W) - 1) % p ** W
    sg = sigma_val(lam, y, p, W - 5)
    # kappa = sg / (L lam^n0)  -> Bewertung, Einheitsteil
    den = L * pow(lam, n0, p ** W) if n0 >= 0 else L * pow(pow(lam // p ** vp(lam, p), -1, p ** W), 1, p ** W)
    if n0 >= 0:
        vden = vp(L, p) + n0 * vp(lam, p)
        uden = (L // p ** vp(L, p)) * pow(lam // p ** vp(lam, p), n0, p ** W)
    else:  # n0 = -1: kappa = lam * sigma(-lam) / L
        vden = vp(L, p) - vp(lam, p)
        uden = (L // p ** vp(L, p)) * pow(lam // p ** vp(lam, p), -1, p ** W)
    vs = vp(sg, p)
    us = sg // p ** vs
    prec = Nw - 5 - vs - 10
    ku = us * pow(uden % p ** prec, -1, p ** prec) % p ** prec
    return ku, vs - vden, prec, n0, [z == 1 for z in zetas[:n0i + 1]]


def kappa_direct(a, p, N2):
    """wie Bericht 04 E2c: kappa_n = (a^^n - alpha)/lam^n bei grosser Hoehe."""
    q = p ** N2
    al = ar_lib.tower_mod(a, 10 * N2, q)
    lam = al * ar_lib.padic_log(a, p, N2) % q
    V = vp(lam, p)
    lu = lam // p ** V
    best = None
    for n in range(2, 400):
        d = (ar_lib.tower_mod(a, n, q) - al) % q
        v = vp(d, p)
        if v + 2 * V + 10 > N2:
            break
        prec = min(N2 - v, v + V) - 2   # Fehler kappa_n - kappa ~ p^(v(d) + V) relativ
        if prec < 1 or (best is not None and prec <= best[2]):
            continue
        ku = (d // p ** v) * pow(pow(lu, n, p ** prec), -1, p ** prec) % p ** prec
        best = (ku, v - n * V, prec)
    return best


bad = 0; cnt = 0
for p in (2, 3, 5, 7, 11, 13):
    for a in range(2, 60):
        if a % p == 0:
            continue
        ku, kv, prec, n0, patt = kappa_closed(a, p, N)
        kd = kappa_direct(a, p, 130)
        P_ = min(prec, kd[2], 25)
        ok = (kv == kd[1]) and (ku - kd[0]) % p ** P_ == 0
        cnt += 1
        if not ok:
            bad += 1
            print("  ABW", p, a, "v:", kv, kd[1], "Uebereinst.:", vp(ku - kd[0], p), "von", P_)
        elif a < 8:
            print(f"  p={p} a={a} n0={n0} v(kappa)={kv}  Einheitsteil mod p^{P_} gleich  (Stellen: geschlossen {prec}, direkt {kd[2]})")
print("Faelle", cnt, "Abweichungen", bad, "Zeit", round(time.time() - t0, 1), "s")
