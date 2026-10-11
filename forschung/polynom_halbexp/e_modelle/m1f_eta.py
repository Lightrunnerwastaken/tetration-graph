"""Modell 1f: erster Newton-Schritt von T_N f_M aus, η_N = |J⁻¹F(T_N f_M)| (Norm Σ|δ_k|ρ^k), für die elliptische Möbius-Abbildung.
Fällt η_N, liegt (heuristisch) eine Lösung nahe T_N f_M; wächst η_N, gibt es dort keine."""
import sys
import mpmath as mp
from flint import arb, ctx
from kern import F_J, g_moebius_ell

mp.mp.dps = 50
CS = [mp.mpf(x) for x in sys.argv[1].split(",")]
Ns = [int(x) for x in sys.argv[2].split(",")]
rho = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf("0.5")
ctx.prec = 1200
for c in CS:
    tau = mp.tan(c / 2)
    zeile = []
    for N in Ns:
        g = [arb(x) for x in g_moebius_ell(c)(N)]
        fk = [tau] + [(1 + tau ** 2) * tau ** (k - 1) for k in range(1, N + 1)]
        a0 = [arb(mp.nstr(x, 300)) for x in fk]
        F, J = F_J(a0, g, N)
        d = J.solve(F)
        eta = mp.fsum(abs(mp.mpf(d[k, 0].mid().str(30, radius=False))) * rho ** k for k in range(N + 1))
        zeile.append(f"N={N}: {mp.nstr(eta, 3)}")
    print(f"c={mp.nstr(c, 4)} tau={mp.nstr(tau, 5)} rho={rho}: " + "  ".join(zeile), flush=True)
