"""Strang E — Kern: Abschneide-Verfahren für f∘f = g (oder f∘f∘f = g) bei beliebigem g.

p_N(x) = Σ_{k≤N} a_k x^k mit [x^k] p∘p = g_k (k = 0..N), Zweig stetig von N = 1 aus (Newton, Start p_{N−1} mit a_N = 0).
Gerechnet mit arb (python-flint). g wird als Funktion gk(N) -> [g_0..g_N] (arb) übergeben.
"""
import sys

import mpmath as mp
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")


def potenzen(f, N):
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    return P


def kompos(a, P, N):
    """Σ a_m P[m] (= p∘q, wenn P die Potenzen von q sind)."""
    s = arb_poly([0])
    for m in range(len(a)):
        s += a[m] * P[m]
    return s.truncate(N + 1)


def abl_kompos(a, P, N):
    """p'∘q = Σ m a_m P[m−1]."""
    s = arb_poly([0])
    for m in range(1, len(a)):
        s += m * a[m] * P[m - 1]
    return s.truncate(N + 1)


def F_J(a, g, N, ordnung=2):
    """Residuum F_k = [x^k] p^{∘ordnung} − g_k und Jacobi-Matrix (arb)."""
    p = arb_poly(a)
    P = potenzen(p, N)
    F, J = arb_mat(N + 1, 1), arb_mat(N + 1, N + 1)
    if ordnung == 2:
        pp = kompos(a, P, N)
        dpp = abl_kompos(a, P, N)
        for k in range(N + 1):
            F[k, 0] = pp[k] - g[k]
            for j in range(N + 1):
                J[k, j] = P[j][k] + (dpp[k - j] if j <= k else 0)
    elif ordnung == 3:
        q = kompos(a, P, N)                     # p∘p
        Q = potenzen(q, N)
        ppp = kompos(a, Q, N)                   # p∘p∘p
        dpq = abl_kompos(a, Q, N)               # p'∘(p∘p)
        dpp = abl_kompos(a, P, N)               # p'∘p
        for k in range(N + 1):
            F[k, 0] = ppp[k] - g[k]
        for j in range(N + 1):
            xj = arb_poly([0] * j + [1])
            spalte = (Q[j] + (dpq * (P[j] + (dpp * xj).truncate(N + 1))).truncate(N + 1))
            for k in range(N + 1):
                J[k, j] = spalte[k]
    else:
        raise ValueError(ordnung)
    return F, J


def newton(a, g, N, ordnung=2, maxit=100, ziel=None):
    """Newton für F(a) = ziel (ziel: arb_mat-Spalte oder None = 0)."""
    alt = None
    for it in range(maxit):
        F, J = F_J(a, g, N, ordnung)
        if ziel is not None:
            F = F - ziel
        d = J.solve(F)
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
        nd = max(abs(d[i, 0].mid()) for i in range(N + 1))
        if nd < arb(2) ** (-ctx.prec + 60):
            return a, it
        # Stagnation auf Rundungsniveau (schlechte Kondition): Schritt winzig und kaum noch kleiner
        if alt is not None and nd < arb(2) ** (-ctx.prec // 2) and nd > alt / 4:
            return a, it
        alt = nd
        if it > 30 and nd > 1e3:
            break
    raise RuntimeError(f"Newton ohne Konvergenz bei N = {N}")


def homotopie(a, g, N, ordnung=2, schritte=8, tief=0):
    """Weg F(a) = (1−t)·F(a_start), t: 0 → 1 (Fortsetzung mit Schrittweitensteuerung)."""
    F0, _ = F_J(a, g, N, ordnung)
    F0 = arb_mat([[arb(F0[i, 0].mid())] for i in range(N + 1)])
    t, h = arb(0), arb(1) / schritte
    while t < 1:
        tn = t + h if t + h < 1 else arb(1)
        try:
            a_neu, _ = newton(a, g, N, ordnung, maxit=40, ziel=(1 - tn) * F0)
            a, t = a_neu, tn
            h = h * 2
        except (RuntimeError, ZeroDivisionError, ValueError):
            h = h / 4
            if h < arb(2) ** -40:
                raise RuntimeError(f"Homotopie bleibt stecken bei N = {N}, t = {t.mid()}")
    return a, -1


def zweig(gk, nmax, bits=500, ordnung=2, a1=None, start=None, ausgabe=None, still=True):
    """{N: [a_0..a_N] (mpmath)} für N = 1..nmax. a1: Wahl von a_1 bei N = 1 (sonst positive Wurzel)."""
    alt, ctx.prec = ctx.prec, bits
    try:
        g = gk(nmax)
        g = [arb(x) for x in g]
        if start is None:
            if ordnung == 2:
                a_1 = arb(a1) if a1 is not None else g[1].sqrt()
                a_0 = g[0] / (1 + a_1)
            else:
                a_1 = arb(a1) if a1 is not None else g[1] ** (arb(1) / 3)
                a_0 = g[0] / (1 + a_1 + a_1 * a_1)
            a = [a_0, a_1]
        else:
            a = [arb(x) for x in start]
        aus = {}
        N0 = len(a) - 1
        for N in range(N0, nmax + 1):
            a = a + [arb(0)] * (N + 1 - len(a))
            try:
                a, it = newton(a, g[:N + 1], N, ordnung)
            except (RuntimeError, ZeroDivisionError, ValueError):
                a, it = homotopie(a, g[:N + 1], N, ordnung)
                if not still:
                    print(f"  N = {N}: Newton direkt gescheitert, Homotopie benutzt", flush=True)
            aus[N] = [mp.mpf(c.mid().str(max(30, int(bits * 0.3011) + 5), radius=False)) for c in a]
            if not still:
                print(N, it, mp.nstr(aus[N][0], 15), flush=True)
            if ausgabe is not None:
                ausgabe(N, aus[N])
        return aus
    finally:
        ctx.prec = alt


def horner(c, z):
    s = mp.mpc(0) if isinstance(z, mp.mpc) else mp.mpf(0)
    for x in reversed(c):
        s = s * z + x
    return s


# ---------- Modell-Funktionen g (Taylor-Koeffizienten als arb) ----------

def g_moebius_ell(c):
    """g(x) = tan(arctan x + c) = (x + t)/(1 − t x), t = tan c."""
    def gk(N):
        t = arb(c).tan() if not isinstance(c, arb) else c.tan()
        return [t] + [(1 + t * t) * t ** (k - 1) for k in range(1, N + 1)]
    return gk


def g_moebius_hyp(c):
    """g(x) = tanh(artanh x + c) = (x + T)/(1 + T x), T = tanh c."""
    def gk(N):
        T = arb(c).tanh() if not isinstance(c, arb) else c.tanh()
        return [T] + [(1 - T * T) * (-T) ** (k - 1) for k in range(1, N + 1)]
    return gk


def g_exp_basis(b):
    """g(x) = b^x."""
    def gk(N):
        lb = arb(b).log()
        out, t = [], arb(1)
        for k in range(N + 1):
            out.append(t)
            t = t * lb / (k + 1)
        return out
    return gk


def g_expm1():
    def gk(N):
        out, t = [arb(0)], arb(1)
        for k in range(1, N + 1):
            out.append(t)
            t = t / (k + 1)
        return out
    return gk
