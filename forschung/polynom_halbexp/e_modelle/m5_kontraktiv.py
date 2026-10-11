"""Modell 5 (Satz K, kontraktiver Fall): g = f∘f mit f(x) = 1/(3 − x) (Möbius, Fixpunkte (3∓√5)/2, anziehend bei 0.382).
‖f‖_{W_r} = Σ|f_k| r^k = 1/(3 − r) < r für 0.382 < r < 2.618 ⇒ Satz K greift: p_N → f geometrisch.
Ausgabe: Fehler a_0, p_N(1), ‖p_N − T_N f‖_{W_r} für r = 1, 2; Newton-Schritt η_N; |A_N^{-1}| in W_1."""
import sys
import mpmath as mp
from flint import arb, ctx, arb_mat
from kern import zweig, horner, F_J

mp.mp.dps = 80
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
f = lambda x: 1 / (3 - x)
g = lambda x: (3 - x) / (8 - 3 * x)
fk = [mp.mpf(1) / 3 ** (k + 1) for k in range(NMAX + 1)]
gk_mp = mp.taylor(g, 0, NMAX)
gk = lambda N: [arb(mp.nstr(x, 75)) for x in gk_mp[:N + 1]]
loes = zweig(gk, NMAX, bits=600)
for N in (1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60, 80):
    if N in loes:
        a = loes[N]
        w = {r: mp.fsum(abs(a[k] - fk[k]) * mp.mpf(r) ** k for k in range(N + 1)) for r in (1, 2)}
        print(f"N={N:3d} a0-f(0)={mp.nstr(a[0] - fk[0], 4):>11} p(1)-f(1)={mp.nstr(horner(a, mp.mpf(1)) - f(mp.mpf(1)), 4):>11}"
              f"  |p-T_Nf|_W1={mp.nstr(w[1], 4):>10}  |p-T_Nf|_W2={mp.nstr(w[2], 4):>10}")
# gleichmässige Schranke für die Inverse der Linearisierung in W_1 (gleiche Gewichte)
ctx.prec = 600
for N in (10, 20, 40, 60):
    if N > NMAX:
        break
    a = [arb(mp.nstr(x, 150)) for x in fk[:N + 1]]
    _, J = F_J(a, [arb(mp.nstr(x, 75)) for x in gk_mp[:N + 1]], N)
    Ji = J.inv()
    for r in (1, 2):
        rr = mp.mpf(r)
        nrm = max(mp.fsum(abs(mp.mpf(Ji[k, j].mid().str(20, radius=False))) * rr ** k for k in range(N + 1)) / rr ** j for j in range(N + 1))
        print(f"N={N}: |J(T_N f)^-1| in W_{r} = {mp.nstr(nrm, 5)}")
