"""k01 — Ist die SVD-freie Majorante G_N(t) = ‖J̃_N^{-1} ẽ_t‖ auf γ = Linse(Γin) ∪ Γout polynomial in N?

Hintergrund (Strang K, Satz K1 in Arbeit): Für jede Eingabe W, analytisch auf einer Umgebung von Int(γ)‾, gilt
  J_N^{-1} T_N W = (1/2πi) ∮_γ W(t) · J_N^{-1} T_N k_t dt,  also ‖J_N^{-1}T_N W‖ ≤ (1/2π) ∮_γ |W(t)| G_N(t) |dt|/|t|.
Mit W = J_∞ k_η folgt eine Pol-Antwort-Schranke ohne Eckformel D6 — gebraucht wird nur G_N ≤ C N^κ auf γ.

Punkte: Γout (|t| = ρ, Winkel 0°…76°), Linse über Γin (Winkel arg L + τ, Radius ρ·exp(Λ(τ)/x0) mit Ws Profil Λ, W §18),
Kontrolle: Γin auf dem Kreis (180°, 120°), L selbst.
Aufruf: python k01_G_kurve.py kubisch 100,200,300,400
"""
import math
import sys
import time

from flint import ctx

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from k_kern import J_tilde, to_acb_mat, cauchy_rhs, spaltennormen, Lc, RHO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")

argL = math.atan2(Lc.imag, Lc.real)
# Ws Kreisprofil Λ(τ) (alle Eckmoden, W §18): τ in Grad ab arg L Richtung Γin
TAU = [0, 2, 4, 8, 15, 30, 60, 103.38]
LAM = [0.0, 0.0009, 0.0025, 0.0049, 0.0069, 0.0091, 0.0109, 0.0116]


def Lam(tau):
    for i in range(len(TAU) - 1):
        if TAU[i] <= tau <= TAU[i + 1]:
            w = (tau - TAU[i]) / (TAU[i + 1] - TAU[i])
            return LAM[i] + w * (LAM[i + 1] - LAM[i])
    return LAM[-1]


def punkte():
    P = []
    for deg in (0, 30, 60, 70, 74, 76):
        P.append(("Gout %d°" % deg, RHO * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))))
    P.append(("L", Lc))
    for x0 in (0.7, 0.95):
        for tau in (2, 4, 8, 15, 30, 60, 103.38):
            s = Lam(tau) / x0
            ph = argL + math.radians(tau)
            P.append(("Linse x0=%.2f τ=%g°" % (x0, tau), RHO * math.exp(s) * complex(math.cos(ph), math.sin(ph))))
    for deg in (120, 180):
        P.append(("Gin-Kreis %d°" % deg, RHO * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))))
    return P


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "kubisch"
    Ns = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "100,200,300").split(",")]
    P = punkte()
    tab = {}
    for N in Ns:
        t0 = time.time()
        J = J_tilde(name, N)
        Ja = to_acb_mat(J)
        B = cauchy_rhs(N, [p[1] for p in P])
        X = Ja.solve(B, algorithm="approx")
        G = spaltennormen(X)
        tab[N] = G
        print("N = %d  (%.1f s, prec %d)" % (N, time.time() - t0, ctx.prec), flush=True)
    print()
    kopf = "%-26s" % "Punkt" + "".join("%12s" % ("ln G/%d" % N) for N in Ns) + "   ln G − ½ln N (letztes)"
    print(kopf)
    for i, (nm, t) in enumerate(P):
        row = "%-26s" % nm + "".join("%12.3f" % math.log(tab[N][i]) for N in Ns)
        row += "   %8.3f" % (math.log(tab[Ns[-1]][i]) - 0.5 * math.log(Ns[-1]))
        print(row)
    if len(Ns) >= 2:
        print("\nRate (ln G(N2) − ln G(N1))/(N2 − N1) zwischen den letzten zwei N:")
        a, b = Ns[-2], Ns[-1]
        for i, (nm, t) in enumerate(P):
            print("%-26s %+.5f" % (nm, (math.log(tab[b][i]) - math.log(tab[a][i])) / (b - a)))


if __name__ == "__main__":
    main()
