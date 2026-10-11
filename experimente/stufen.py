"""Höhere Stufen: Taylorreihen der Kneser-Tetration aus Kreiswerten, Fixpunkte von T_b, Poincaré-Funktion an einem Fixpunkt von T_b.

Alle Kneser-Werte über gemeinsam.kneser (fatou.gp-Fork, Zwischenspeicher). Polynome als Koeffizientenlisten [a_0, a_1, …] in t = x − x0.
"""
import mpmath as mp

from gemeinsam import kneser


def basis_text(b):
    """Basis als GP-Ausdruck; 'exp(1)' bleibt, Zahlen mit 45 Stellen."""
    return b if isinstance(b, str) else mp.nstr(mp.mpf(b), 45)


def taylor(basis, x0, r, N=96, K=40, art="sexp"):
    """Taylor-Koeffizienten a_0..a_K von T_b um x0 aus N Werten auf |z − x0| = r (Cauchy, diskret)."""
    zs = [x0 + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
    w = kneser(art, basis_text(basis), zs)
    a = []
    for m in range(K + 1):
        s = mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N
        a.append(mp.re(s) / r ** m)
    return a


def p_wert(a, t):
    s = mp.mpf(0)
    for c in reversed(a):
        s = s * t + c
    return s


def p_abl(a):
    return [k * a[k] for k in range(1, len(a))]


def p_verschieben(a, t0):
    """Koeffizienten um x0 + t0 statt x0."""
    K = len(a) - 1
    return [mp.fsum(a[m] * mp.binomial(m, k) * t0 ** (m - k) for m in range(k, K + 1)) for k in range(K + 1)]


def p_nullstelle(a, t=0, schritte=60):
    """Newton für p(t) = 0."""
    da = p_abl(a)
    for _ in range(schritte):
        d = p_wert(a, t) / p_wert(da, t)
        t -= d
        if abs(d) < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
    return t


def fixpunkt(basis, x0, r, **kw):
    """Fixpunkt L von T_b nahe x0 und Taylorreihe um L (a[1] = Multiplikator)."""
    a = taylor(basis, x0, r, **kw)
    g = list(a)                                # p(t) − (x0 + t) = (a_0 − x0) + (a_1 − 1) t + a_2 t² + …
    g[0], g[1] = a[0] - x0, a[1] - 1
    t = p_nullstelle(g)
    return x0 + t, p_verschieben(a, t)


def poincare(a, K=30):
    """Poincaré-Funktion Ψ an L für f mit Taylorreihe a um L (a_0 = L, a_1 = λ): Ψ(λw) = f(Ψ(w)), Ψ(0) = L, Ψ'(0) = 1.

    Gibt [ψ_1 = 1, ψ_2, …, ψ_K] (Koeffizienten von Ψ − L) zurück.
    """
    lam = a[1]
    psi = [mp.mpf(0), mp.mpf(1)]
    for n in range(2, K + 1):
        # [w^n] Σ_{k≥2} a_k Q^k mit Q = Σ_{j<n} ψ_j w^j (ψ_n trägt dazu nicht bei)
        Q = psi + [mp.mpf(0)] * (n + 1 - len(psi))
        potenz = list(Q)
        s = mp.mpf(0)
        for k in range(2, n + 1):
            neu = [mp.mpf(0)] * (n + 1)
            for i, pi in enumerate(potenz):
                if pi == 0:
                    continue
                for j in range(1, n + 1 - i):
                    neu[i + j] += pi * Q[j]
            potenz = neu
            if k < len(a):
                s += a[k] * potenz[n]
        psi.append(s / (lam ** n - lam))
    return psi[1:]


def psi_wert(L, psi, w):
    s = mp.mpf(0)
    for c in reversed(psi):
        s = (s + c) * w
    return L + s


def pentation(basis, x0, r, M=8, K=30, n_slog=7, N=128):
    """Reguläre Pentation zur Basis b am Fixpunkt L_b von T_b nahe x0: (L, λ, ψ, c, slog-Kette, P).

    P(zs, M) rechnet P_b(z) = T_b^M(Ψ(c λ^{z−M})) für einen Stapel (reell oder komplex).
    """
    L, a = fixpunkt(basis, x0, r, N=N, K=40)
    lam = a[1]
    psi = poincare(a, K=K)
    kette = [mp.mpf(1)]                                   # 1, slog(1) = 0, slog(0) = −1, …
    for _ in range(n_slog):
        kette.append(mp.re(kneser("slog", basis_text(basis), [kette[-1]])[0]))
    c = lam ** n_slog * psi_umkehr(L, psi, kette[-1])

    def P(zs, M=M):
        reell = all(mp.im(mp.mpc(z)) == 0 for z in zs)
        ys = [psi_wert(L, psi, c * lam ** (z - M)) for z in zs]
        for _ in range(M):
            ys = kneser("sexp", basis_text(basis), ys)
            if reell:
                ys = [mp.re(v) for v in ys]
        return ys
    return L, lam, psi, c, kette, P


def psi_umkehr(L, psi, z, w=None):
    """w mit Ψ(w) = z (z nahe L)."""
    w = z - L if w is None else w
    dpsi = [k * c for k, c in enumerate(psi, start=1)]
    for _ in range(60):
        f = psi_wert(L, psi, w) - z
        d = mp.fsum(c * w ** (k - 1) for k, c in enumerate(dpsi, start=1))
        dw = f / d
        w -= dw
        if abs(dw) < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
    return w
