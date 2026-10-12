"""k02 — Ist die Ecke von G_N (Wert in L und nahe L auf γ) lokal bestimmt?

Vergleich von Modellen mit denselben Fixpunktdaten (L, μ) und verschiedener globaler Gestalt (Gs jet_ε: λ ändert sich um ≈ 0.0097),
dazu das radiale Quadrat (gleiches L, gleiches |μ|, arg μ = 31.4° statt 38.3°).
Ist G_N(L)/√N (und die Profile nahe L) für jet_±0.01 gleich wie kubisch, ist die Eck-Amplitude in der G-Form lokal (Gegenprobe zu W §27.5).
Aufruf: python k02_ecke_lokal.py 300,500,700
"""
import math
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from k_kern import J_tilde, to_acb_mat, cauchy_rhs, spaltennormen, Lc, RHO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
argL = math.atan2(Lc.imag, Lc.real)


def punkte():
    P = [("L", Lc)]
    for deg in (76.0, 75.0, 70.0, 0.0):
        P.append(("Gout %g°" % deg, RHO * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))))
    for tau, s in ((0.5, 0.0004), (1.0, 0.0007), (2.0, 0.0013), (15.0, 0.0099)):   # Linse x0 = 0.7 (Ws Profil)
        ph = argL + math.radians(tau)
        P.append(("Linse τ=%g°" % tau, RHO * math.exp(s) * complex(math.cos(ph), math.sin(ph))))
    # Strahl nach aussen (Schlitzrichtung) in Abstand u·ρ
    for u in (0.002, 0.01):
        P.append(("Schlitz u=%g" % u, Lc * (1 + u)))
    return P


def main():
    Ns = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "300,500").split(",")]
    modelle = sys.argv[2].split(",") if len(sys.argv) > 2 else ["kubisch", "jet0.01", "jet-0.01", "radial"]
    P = punkte()
    res = {}
    for name in modelle:
        for N in Ns:
            t0 = time.time()
            J = to_acb_mat(J_tilde(name, N))
            X = J.solve(cauchy_rhs(N, [p[1] for p in P]), algorithm="approx")
            res[(name, N)] = spaltennormen(X)
            print("%s N=%d (%.1f s)" % (name, N, time.time() - t0), flush=True)
    print("\nln G − ½ ln N:")
    for i, (nm, t) in enumerate(P):
        print("%-16s" % nm + "  ".join(
            "%s:" % name[:8] + " ".join("%7.3f" % (math.log(res[(name, N)][i]) - 0.5 * math.log(N)) for N in Ns)
            for name in modelle))


if __name__ == "__main__":
    main()
