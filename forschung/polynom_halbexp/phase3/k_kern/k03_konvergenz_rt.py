"""k03 — Konvergiert die Abschnittsantwort r_t^{(N)} = J_N^{-1} T_N k_t koeffizientenweise (gegen s_t = J_∞^{-1} k_t), und wie schnell?

Für t auf γ (Γout, Linse, L) wird r_t^{(N)} mit T_N r_t^{(N')} (N' gross, Ersatz für s_t) verglichen:
  δ_N(t) := ‖r_t^{(N)} − T_N r_t^{(N')}‖_{ℓ²_ρ},  dazu ‖T_{>N} r_t^{(N')}‖ (Schwanz von s_t) und die Koeffizientenprofile.
Aufruf: python k03_konvergenz_rt.py kubisch 100,200,300,400 800
"""
import math
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, ctx  # noqa: E402
from k_kern import J_tilde, to_acb_mat, cauchy_rhs, Lc, RHO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
argL = math.atan2(Lc.imag, Lc.real)


def punkte():
    P = []
    for deg in (0.0, 60.0, 74.0):
        P.append(("Gout %g°" % deg, RHO * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))))
    P.append(("L", Lc))
    for tau, s in ((2.0, 0.0013), (15.0, 0.0099), (103.38, 0.0166)):
        ph = argL + math.radians(tau)
        P.append(("Linse τ=%g°" % tau, RHO * math.exp(s) * complex(math.cos(ph), math.sin(ph))))
    return P


def loese(name, N, P):
    J = to_acb_mat(J_tilde(name, N, bits=int(1.0 * N) + 300))
    return J.solve(cauchy_rhs(N, [p[1] for p in P]), algorithm="approx")


def main():
    name = sys.argv[1]
    Ns = [int(x) for x in sys.argv[2].split(",")]
    Nref = int(sys.argv[3])
    P = punkte()
    t0 = time.time()
    Xref = loese(name, Nref, P)
    print("Referenz N' = %d (%.1f s)" % (Nref, time.time() - t0), flush=True)
    for c, (nm, t) in enumerate(P):
        # Schwanzprofil der Referenzlösung: max |coef| in Fenstern
        prof = []
        for a, b in ((0, 50), (100, 150), (200, 250), (300, 350), (400, 450), (600, 650)):
            if b <= Nref:
                m = max(float(abs(Xref[k, c]).mid()) for k in range(a, b))
                prof.append("%d:%.1e" % (a, m))
        print("%-14s Ref-Koeff.-Profil (max|.| je Fenster, ρ-skaliert): %s" % (nm, " ".join(prof)))
    print()
    for N in Ns:
        X = loese(name, N, P)
        zeile = []
        for c, (nm, t) in enumerate(P):
            d = arb(0)
            for k in range(N + 1):
                d += abs(X[k, c] - Xref[k, c]) ** 2
            tail = arb(0)
            for k in range(N + 1, Nref + 1):
                tail += abs(Xref[k, c]) ** 2
            zeile.append((nm, float(d.sqrt().mid()), float(tail.sqrt().mid())))
        print("N = %d:" % N)
        for nm, d, tl in zeile:
            print("   %-14s ‖r^(N) − T_N r^(N')‖ = %.3e    ‖Schwanz_{>N} r^(N')‖ = %.3e" % (nm, d, tl))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
