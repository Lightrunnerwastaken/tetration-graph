"""X03 — Struktur von J_N für Möbius-Modelle: Kompositionsteil = (untere Dreiecks-Toeplitz) × (SL(2)-Darstellung π_N(g)).

Behauptung (Rechnung): für f(x) = (ax+b)/(cx+d) und P ∈ P_N gilt P∘f = (cx+d)^{−N}·(π_N(g)P)(x) mit
(π_N(g)P)(x) := (cx+d)^N P((ax+b)/(cx+d)) ∈ P_N (irreduzible (N+1)-dim. Darstellung von GL(2), Matrixelemente = Krawtchouk-Polynome).
Also T_N[P∘f] = L_N(r)·π_N(g)P, r = (cx+d)^{−N}, L_N(r) = untere Dreiecks-Toeplitz-Matrix der Taylor-Koeffizienten von r bei 0.
π_N(g) ist diagonal in der Fixpunktbasis e_j = (x−L₁)^j (x−L₂)^{N−j} mit Eigenwerten (a−cL₁)^j (a−cL₂)^{N−j}.
Prüfung: zufällige komplexe Möbius-Abbildung, N = 8, 20 (mpmath, 50 Stellen).
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 50


def taylor_komp(P, a, b, c, d, N):
    """T_N[P∘f] direkt: Taylor-Reihe von f bei 0, Potenzen, abgeschnitten."""
    # f = (ax+b)/(cx+d) = Σ f_k x^k
    fk = [mp.mpc(0)] * (N + 1)
    # 1/(cx+d) = (1/d) Σ (−c/d)^k x^k
    inv = [(-c / d) ** k / d for k in range(N + 1)]
    for k in range(N + 1):
        fk[k] = a * (inv[k - 1] if k >= 1 else 0) + b * inv[k]
    out = [mp.mpc(0)] * (N + 1)
    pot = [mp.mpc(1)] + [mp.mpc(0)] * N
    for j in range(len(P)):
        for k in range(N + 1):
            out[k] += P[j] * pot[k]
        neu = [mp.mpc(0)] * (N + 1)
        for i in range(N + 1):
            if pot[i] == 0:
                continue
            for m in range(N + 1 - i):
                neu[i + m] += pot[i] * fk[m]
        pot = neu
    return out


def polymul(p, q):
    r = [mp.mpc(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def pi_N(P, a, b, c, d, N):
    """(cx+d)^N P((ax+b)/(cx+d)) = Σ p_j (ax+b)^j (cx+d)^{N−j}."""
    out = [mp.mpc(0)] * (N + 1)
    for j in range(N + 1):
        t = [mp.mpc(1)]
        for _ in range(j):
            t = polymul(t, [b, a])
        for _ in range(N - j):
            t = polymul(t, [d, c])
        for k in range(N + 1):
            out[k] += P[j] * t[k]
    return out


mp.mp.dps = 50
import random
random.seed(3)
for N in (8, 20):
    a, b, c, d = [mp.mpc(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(4)]
    d += 2
    P = [mp.mpc(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(N + 1)]
    links = taylor_komp(P, a, b, c, d, N)
    Q = pi_N(P, a, b, c, d, N)
    r = [d ** (-N) * mp.binomial(-N, k) * (c / d) ** k for k in range(N + 1)]
    rechts = [sum(r[k - j] * Q[j] for j in range(k + 1)) for k in range(N + 1)]
    err = max(abs(x - y) for x, y in zip(links, rechts))
    # Eigenbasis
    disc = mp.sqrt((d - a) ** 2 + 4 * b * c)
    L1, L2 = (a - d + disc) / (2 * c), (a - d - disc) / (2 * c)
    ej = [mp.mpc(1)]
    j = 3
    for _ in range(j):
        ej = polymul(ej, [-L1, 1])
    for _ in range(N - j):
        ej = polymul(ej, [-L2, 1])
    Pe = pi_N(ej, a, b, c, d, N)
    lam = (a - c * L1) ** j * (a - c * L2) ** (N - j)
    err2 = max(abs(x - lam * y) for x, y in zip(Pe, ej))
    print(f"N={N}: |T_N[P∘f] − L_N(r)π_N(g)P| = {mp.nstr(err, 3)};  Eigenvektor e_3: Fehler {mp.nstr(err2, 3)}")
