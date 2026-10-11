"""Strang P, Prüfung 7: unabhängige Stichprobe der numerischen Stütze von (H₀) (B b25): dist(ζ,Ḡ')·max_{|z|=ρ}|Π_N k_ζ(z)|
an der Ecke (Winkelhalbierende) und links bei −ρ, N = 60, 100; eigene Matrix J_N(f_K) aus den 4000 Kneser-Koeffizienten."""
import sys
import time

import mpmath as mp

from reihen import add, compose, deriv, horner, inv, kneser_koeff, mul, shift, fixpunkt_L

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
A = kneser_koeff(4000)
L = fixpunkt_L()
rho, a0 = abs(L), A[0]
R0 = abs(L - a0)
psi = mp.radians(mp.mpf("87.151205"))


def _dist(z):
    # Abstand zu D̄ ∪ D̄(a0,R0) für z ausserhalb: min über beide Kreisscheiben, aber Punkte der Kreise, die in der anderen
    # Scheibe liegen, zählen nicht zum Rand -> Ecke: Abstand zu L/L̄ falls die Projektion im Inneren der anderen Scheibe liegt.
    cands = []
    for (c, r) in [(mp.mpf(0), rho), (a0, R0)]:
        p = c + r * (z - c) / abs(z - c)
        other_c, other_r = (a0, R0) if c == 0 else (mp.mpf(0), rho)
        if abs(p - other_c) >= other_r:
            cands.append(abs(z - p))
    cands += [abs(z - L), abs(z - mp.conj(L))]
    return min(cands)


for N in [60, 100]:
    t0 = time.time()
    Fz = A[:N + 2]
    Fc = shift(A, a0, N + 2)
    m = compose(deriv(Fc) + [0], Fz, N)
    cols = []
    pw = [mp.mpf(1)] + [mp.mpf(0)] * N
    for j in range(N + 1):
        cols.append(add(pw, [0] * j + m[:N + 1 - j], N))
        pw = mul(pw, Fz, N)
    J = mp.matrix([[cols[j][k] for j in range(N + 1)] for k in range(N + 1)])
    zs = [rho * mp.expj(2 * mp.pi * k / 720) for k in range(720)]
    print(f"N = {N} (Matrix {time.time() - t0:.0f} s)")
    for name, basis, richt in [("Ecke L, Winkelhalbierende", L, mp.expj(psi)), ("links −ρ", -rho, mp.mpf(-1))]:
        out = []
        for t in [mp.mpf("0.5"), mp.mpf(2), mp.mpf(8)]:
            zeta = basis + (t / N) * richt
            rhs = add(inv(add([zeta], Fz[:N + 1], N, -1), N), mul(m, [1 / zeta ** (k + 1) for k in range(N + 1)], N), N)
            x = mp.lu_solve(J, mp.matrix(rhs))
            c = [x[k] for k in range(N + 1)]
            mx = max(abs(horner(c, z)) for z in zs)
            out.append(_dist(zeta) * mx)
        print(f"   {name:28s} t = 0.5/2/8: " + " / ".join(mp.nstr(v, 4) for v in out))
    print(f"   ({time.time() - t0:.0f} s)")
