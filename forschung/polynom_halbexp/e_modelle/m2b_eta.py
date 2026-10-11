"""Modell 2b: η_N = |J⁻¹F(T_N f_M)| für die elliptische Möbius-Abbildung mit Fixpunkten A ± iB und Drehwinkel c,
dazu die vorhergesagte Rate ρ = χ(0, f(0)) / χ(g(0), ∞) (chordaler Abstand zur Drehung um A ± iB)."""
import sys
import mpmath as mp
from flint import arb, ctx
from kern import F_J

mp.mp.dps = 60
A, B = mp.mpf(sys.argv[1]), mp.mpf(sys.argv[2])
CS = [mp.mpf(x) for x in sys.argv[3].split(",")]
Ns = [int(x) for x in sys.argv[4].split(",")]
rho = mp.mpf("0.5")
ctx.prec = 1200
for c in CS:
    g = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c)
    fM = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c / 2)
    th0 = mp.atan(-A / B)
    vorh = mp.sin(c / 2) / abs(mp.cos(th0 + c))
    gk = mp.taylor(g, 0, max(Ns))
    fk = mp.taylor(fM, 0, max(Ns))
    zeile = []
    alt = None
    for N in Ns:
        ga = [arb(mp.nstr(x, 58)) for x in gk[:N + 1]]
        a0 = [arb(mp.nstr(x, 58)) for x in fk[:N + 1]]
        F, J = F_J(a0, ga, N)
        d = J.solve(F)
        eta = mp.fsum(abs(mp.mpf(d[k, 0].mid().str(30, radius=False))) * rho ** k for k in range(N + 1))
        zeile.append(f"{N}: {mp.nstr(eta, 3)}")
    print(f"c={mp.nstr(c, 5)} vorhergesagte Rate {mp.nstr(vorh, 5)} | " + "  ".join(zeile), flush=True)
