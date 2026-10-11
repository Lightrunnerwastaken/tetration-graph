"""M08 — Prüfung der exakten Zwei-Punkt-Form der adjungierten Gleichung (beliebiges Polynom f mit f(0) = a₀, f'(0) ≠ 0):
Für v ∈ C^{N+1}, W(t) := Σ_i v_{N−i} t^i, h = f⁻¹ (Zweig mit h(a₀) = 0), Φ_W(y) := W(h(y)) h'(y) (y/h(y))^{N+1}:
    (J^T v)_{N−l} = [t^l]( m·W − PP_{a₀}[Φ_W] ),  l = 0..N,
gleichwertig  (J^T v)_{N−l} = [t^l] { (t−a₀)^{−N−1} ( T⁽⁰⁾_N[(t−a₀)^{N+1} m W] − T^{(a₀)}_N[W(h) h' Λ^{N+1}] ) },  Λ = t(t−a₀)/h(t).
Prüfung mit Zufallsvektor in mpmath (hohe Präzision), radiales Quadrat-Modell und kubisches Modell.
"""
import sys

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80


def reihe_mul(a, b, n):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1) if i < len(a) and k - i < len(b)) for k in range(n)]


def reihe_pot(a, p, n):
    out = [mp.mpf(1)] + [mp.mpf(0)] * (n - 1)
    for _ in range(p):
        out = reihe_mul(out, a, n)
    return out


def jmat_T_v(coef, N, v):
    """(J^T v)_k = Σ_j v_j ([x^j] f^k + m_{j−k}), m = f'∘f."""
    n = N + 1
    f = [mp.mpf(c) for c in coef] + [mp.mpf(0)] * n
    f = f[:n]
    fp = [(i + 1) * f[i + 1] for i in range(len(coef) - 1)]
    P = [[mp.mpf(1)] + [mp.mpf(0)] * (n - 1)]
    for _ in range(N):
        P.append(reihe_mul(P[-1], f, n))
    m = [mp.mpf(0)] * n
    for i, c in enumerate(fp):
        m = [m[k] + c * P[i][k] for k in range(n)]
    out = []
    for k in range(n):
        s = mp.fsum(v[j] * P[k][j] for j in range(n))
        s += mp.fsum(v[j] * m[j - k] for j in range(k, n))
        out.append(s)
    return out, m


def h_reihe_bei_a0(coef, n):
    """Taylor-Koeffizienten von h = f⁻¹ bei y = a₀ (in Potenzen von (y − a₀)), h(a₀) = 0 — Reihenumkehr."""
    g = [mp.mpf(c) for c in coef] + [mp.mpf(0)] * n   # f(x) − a₀ = Σ_{i≥1} g_i x^i
    g = g[:n]
    g[0] = mp.mpf(0)
    # Lagrange-Umkehr per Newton auf Reihen: h(u) mit f(h(u)) − a₀ = u
    h = [mp.mpf(0), 1 / g[1]] + [mp.mpf(0)] * (n - 2)
    for _ in range(int(np.log2(n)) + 3):
        # F(h) = Σ g_i h^i − u, F'(h) = Σ i g_i h^{i−1}
        hp = [mp.mpf(1)] + [mp.mpf(0)] * (n - 1)
        F = [mp.mpf(0)] * n
        Fd = [mp.mpf(0)] * n
        for i in range(1, len(coef)):
            hp_prev = hp
            hp = reihe_mul(hp, h, n)
            F = [F[k] + g[i] * hp[k] for k in range(n)]
            Fd = [Fd[k] + i * g[i] * hp_prev[k] for k in range(n)]
        F[1] -= 1
        # Newton: h ← h − F/F'
        inv = [mp.mpf(0)] * n
        inv[0] = 1 / Fd[0]
        for k in range(1, n):
            inv[k] = -mp.fsum(Fd[i] * inv[k - i] for i in range(1, k + 1)) / Fd[0]
        corr = reihe_mul(F, inv, n)
        h = [h[k] - corr[k] for k in range(n)]
    return h


def pruefe(coef, N, seed=1):
    rng = np.random.default_rng(seed)
    v = [mp.mpf(float(x)) for x in rng.standard_normal(N + 1)]
    lhs, m = jmat_T_v(coef, N, v)
    a0 = mp.mpf(coef[0])
    M = 2 * N + 4
    # Φ_W bei a₀ in u = y − a₀: W(h(u)) h'(u) ((a₀+u)/h(u))^{N+1}
    h = h_reihe_bei_a0(coef, M + 2)
    W = [v[N - i] for i in range(N + 1)]
    # W(h(u))
    Wh = [mp.mpf(0)] * M
    hp = [mp.mpf(1)] + [mp.mpf(0)] * (M - 1)
    for i in range(N + 1):
        Wh = [Wh[k] + W[i] * hp[k] for k in range(M)]
        hp = reihe_mul(hp, h[:M], M)
    hd = [(k + 1) * h[k + 1] for k in range(M)]
    # h(u)/u = h1 + h2 u + ...  →  (a₀+u)^{N+1} (h(u)/u)^{−N−1} u^{−N−1}
    q = h[1:M + 1]
    qinv = [mp.mpf(0)] * M
    qinv[0] = 1 / q[0]
    for k in range(1, M):
        qinv[k] = -mp.fsum(q[i] * qinv[k - i] for i in range(1, k + 1)) / q[0]
    qn = reihe_pot(qinv, N + 1, M)
    yn = [mp.binomial(N + 1, k) * a0 ** (N + 1 - k) if k <= N + 1 else mp.mpf(0) for k in range(M)]
    core = reihe_mul(reihe_mul(reihe_mul(Wh, hd, M), qn, M), yn, M)   # Φ·u^{N+1}, Koeffizienten in u
    # Hauptteil: Σ_{j=1}^{N+1} core[N+1−j] u^{−j};  Taylor bei t = 0: u = t − a₀
    # [t^l] (t − a₀)^{−j} = binom(−j, l) (−a₀)^{−j−l}
    rhs = []
    for l in range(N + 1):
        pp = mp.fsum(core[N + 1 - j] * mp.binomial(-j, l) * (-a0) ** (-j - l) for j in range(1, N + 2))
        mw = mp.fsum(m[i] * W[l - i] for i in range(l + 1))
        rhs.append(mw - pp)
    fehler = max(abs(lhs[N - l] - rhs[l]) for l in range(N + 1))
    groesse = max(abs(x) for x in lhs)
    print(f"N={N}, f={coef}: max|J^T v − Formel| = {mp.nstr(fehler, 5)} (|J^T v| bis {mp.nstr(groesse, 5)})")


th = np.deg2rad(76.618)
eps = 0.31382
pruefe([eps, 1 - 2 * eps * np.cos(th), eps], 12)
pruefe([eps, 1 - 2 * eps * np.cos(th), eps], 25)
pruefe([0.49999214593419417, 0.8739061716168851, 0.2503912665915572, 0.022377341868459022], 15)
