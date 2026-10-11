"""O02 — Elliptische Möbius-Modelle: ist (I + K_N)^{-1} EXAKT eine Kontraktion (‖Y_N‖ ≤ 1)?
Allgemeine Möbius-Abbildung f(x) = (a x + b)/(c x + d) mit Fixpunkten auf |x| = ρ.
Prüfung hochgenau (mpmath, kleine N): Eigenwerte von A*A − I (A = I + K_N in ℓ²_ρ), Rang des Kerns.
Aufruf: python o02_moebius_exakt.py
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60


def taylor_rat(a, b, c, d, n):
    geo = [(-c / d) ** k / d for k in range(n + 1)]
    return [b * geo[0]] + [b * geo[k] + a * geo[k - 1] for k in range(1, n + 1)]


def polmul(p, q, n):
    r = [mp.mpc(0)] * (n + 1)
    for i, x in enumerate(p[:n + 1]):
        if x == 0:
            continue
        for j, y in enumerate(q[:n + 1 - i]):
            r[i + j] += x * y
    return r


def A_matrix(a, b, c, d, rho, N):
    """A = I + K_N, umskaliert (Basis x^k/ρ^k)."""
    fk = taylor_rat(a, b, c, d, N)
    det = a * d - b * c
    # m = det/(c f + d)^2 als Reihe
    q = [d + c * fk[0]] + [c * x for x in fk[1:]]
    qi = [1 / q[0]] + [mp.mpc(0)] * N
    for k in range(1, N + 1):
        s = sum(q[j] * qi[k - j] for j in range(1, k + 1))
        qi[k] = -s / q[0]
    m = [det * x for x in polmul(qi, qi, N)]
    pw = [[mp.mpc(1)] + [mp.mpc(0)] * N]
    for _ in range(N):
        pw.append(polmul(pw[-1], fk, N))
    C = mp.matrix(N + 1, N + 1)
    Mm = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(N + 1):
            s = mp.mpf(rho) ** (k - j)
            C[k, j] = pw[j][k] * s
            Mm[k, j] = (m[k - j] if k >= j else 0) * s
    Kmat = mp.inverse(Mm) * C
    return mp.eye(N + 1) + Kmat


def pruefe(name, a, b, c, d, rho, Ns):
    for N in Ns:
        A = A_matrix(a, b, c, d, rho, N)
        G = A.H * A - mp.eye(N + 1)
        ev = mp.eighe(G, eigvals_only=True) if hasattr(mp, "eighe") else None
        if ev is None:
            ev, _ = mp.eighe(G)
        ev = sorted([mp.re(e) for e in ev])
        nnull = sum(1 for e in ev if abs(e) < mp.mpf(10) ** (-40))
        print(f"{name} N={N}: min eig(A*A − I) = {mp.nstr(ev[0], 5)}, nächste {mp.nstr(ev[1], 5)}, {mp.nstr(ev[2],5)}; "
              f"#|eig| < 1e-40: {nnull}; max eig = {mp.nstr(ev[-1], 5)}", flush=True)


if __name__ == "__main__":
    # symmetrisch: Drehung um ±i, f(x) = (x + τ)/(1 − τ x), ρ = 1
    for cc in (0.6, 1.0, 1.3):
        tau = mp.tan(mp.mpf(cc) / 2)
        pruefe(f"sym c={cc}", mp.mpc(1), tau, -tau, mp.mpc(1), 1, [4, 8, 14])
    # Möbius-exp (E): Drehung um L, L̄ mit g(0) = 1
    L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237",
               "1.33723570143068940890116214319371061253950213846051241687")
    A0, B0 = L.real, L.imag
    cang = mp.atan((1 - A0) / B0) - mp.atan(-A0 / B0)
    h = mp.tan(cang / 2)
    pruefe("moeb_ell", B0 - A0 * h, h * (A0 ** 2 + B0 ** 2), -h, B0 + h * A0, abs(L), [4, 8, 14])
    # Loxodromisch (μ = √L) zum Vergleich
    mu = mp.sqrt(L)
    Lb = mp.conj(L)
    pruefe("moeb_lox", L - Lb * mu, -L * Lb + Lb * mu * L, 1 - mu, -Lb + mu * L, abs(L), [4, 8, 14])
