"""Modell 1g: Stabilität |J_N(T_N f_M)^{-1}| und Konsistenz |F_N(T_N f_M)| in l¹ mit gleichen Gewichten R^k(1+k)^s
für die elliptische Möbius-Drehung (Fixpunkte ±i, Radius 1). Unter/über der Schwelle c = π/3."""
import sys
import mpmath as mp
from flint import arb, ctx
from kern import F_J, g_moebius_ell

mp.mp.dps = 40
CS = [mp.mpf(x) for x in sys.argv[1].split(",")]
Ns = [int(x) for x in sys.argv[2].split(",")]
RS = [mp.mpf(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [mp.mpf(1)]
S = int(sys.argv[4]) if len(sys.argv) > 4 else 0
ctx.prec = 1000
for c in CS:
    tau = mp.tan(c / 2)
    for R in RS:
        w = lambda k: R ** k * (1 + k) ** S
        zeile = []
        for N in Ns:
            g = [arb(x) for x in g_moebius_ell(c)(N)]
            fk = [tau] + [(1 + tau ** 2) * tau ** (k - 1) for k in range(1, N + 1)]
            a0 = [arb(mp.nstr(x, 250)) for x in fk]
            F, J = F_J(a0, g, N)
            Ji = J.inv()
            Jm = [[mp.mpf(Ji[k, j].mid().str(25, radius=False)) for j in range(N + 1)] for k in range(N + 1)]
            nrm = max(mp.fsum(abs(Jm[k][j]) * w(k) for k in range(N + 1)) / w(j) for j in range(N + 1))
            res = mp.fsum(abs(mp.mpf(F[k, 0].mid().str(25, radius=False))) * w(k) for k in range(N + 1))
            d = J.solve(F)
            eta = mp.fsum(abs(mp.mpf(d[k, 0].mid().str(25, radius=False))) * w(k) for k in range(N + 1))
            zeile.append(f"N={N}: |J⁻¹|={mp.nstr(nrm, 4)} |F|={mp.nstr(res, 3)} η={mp.nstr(eta, 3)}")
        print(f"c={mp.nstr(c, 4)} R={R} s={S}: " + "  ".join(zeile), flush=True)
