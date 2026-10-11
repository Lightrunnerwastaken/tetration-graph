"""Modell 1c: Gibt es für jedes N eine Lösung nahe dem Taylor-Polynom T_N f_M der Möbius-Wurzel (Newton direkt von T_N f_M)?
Dazu: Singulärwerte der gewichteten Jacobi-Matrix in T_N f_M und Anteil des Residuums in den fast-singulären Richtungen."""
import sys
import mpmath as mp
from flint import arb, ctx
from kern import F_J, newton, g_moebius_ell, horner

mp.mp.dps = 60
c = mp.mpf(sys.argv[1])
Ns = [int(x) for x in sys.argv[2].split(",")]
ctx.prec = 800
tau = mp.tan(c / 2)
fM = lambda z: (z + tau) / (1 - tau * z)
gk = g_moebius_ell(c)
print(f"c = {c}, tau = {mp.nstr(tau, 6)}")
for N in Ns:
    g = [arb(x) for x in gk(N)]
    fk = [tau] + [(1 + tau ** 2) * tau ** (k - 1) for k in range(1, N + 1)]
    a0 = [arb(mp.nstr(x, 200)) for x in fk]
    F, J = F_J(a0, g, N)
    res = max(abs(F[k, 0].mid()) for k in range(N + 1))
    try:
        a, it = newton(a0, g, N)
        am = [mp.mpf(x.mid().str(80, radius=False)) for x in a]
        d = horner(am, mp.mpf("0.5")) - fM(mp.mpf("0.5"))
        dk = max(abs(am[k] - fk[k]) for k in range(N + 1))
        print(f"N={N:3d} |F(T_N f)|_max={mp.nstr(mp.mpf(res.str(20, radius=False)), 3)}  Newton {it} Schritte: p(.5)-f(.5)={mp.nstr(d, 3)}  max|da_k|={mp.nstr(dk, 3)}", flush=True)
    except Exception as e:
        print(f"N={N:3d} Newton von T_N f_M gescheitert: {e}")
