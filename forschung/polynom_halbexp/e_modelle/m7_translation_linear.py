"""Modell 7 (Satz T): lineares Translationsmodell u(x + s) + u(x) = h(x), abgeschnitten: u_N ∈ P_N mit T_N[u_N(·+s) + u_N] = T_N h.
Exakt: u_N = φ_s(D) T_N h, φ_s(z) = 1/(1 + e^{sz}). Test mit h = e^{λx}: Grenzwert e^{λx}/(1 + e^{λs}) genau für |λ| < π/s.
(Das ist die Linearisierung der halben Iterierten von g = x + 2s + εh an f = x + s.)"""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")

mp.mp.dps = 50
s = mp.mpf("0.5")                      # π/s = 2π
for lam in (mp.mpf(3), mp.mpf(6), mp.mpf("6.2"), mp.mpf("6.4"), mp.mpf(8), mp.mpc(0, 5), mp.mpc(0, 7)):
    zeile = []
    for N in (10, 20, 40, 80, 160):
        # Matrix von (S + I) auf P_N in Monombasis, rechte Seite T_N e^{λx}
        M = mp.matrix(N + 1, N + 1)
        for j in range(N + 1):            # Spalte j: (x+s)^j + x^j
            for i in range(j + 1):
                M[i, j] += mp.binomial(j, i) * s ** (j - i)
            M[j, j] += 1
        rhs = mp.matrix([lam ** k / mp.factorial(k) for k in range(N + 1)])
        u = mp.lu_solve(M, rhs)
        wert = u[0]
        soll = 1 / (1 + mp.exp(lam * s))
        zeile.append(f"N={N}: {mp.nstr(abs(wert - soll), 3)}")
    print(f"λ = {mp.nstr(lam, 4)} (|λ|/(π/s) = {mp.nstr(abs(lam) * s / mp.pi, 3)}): |u_N(0) − u(0)|  " + "  ".join(zeile))
