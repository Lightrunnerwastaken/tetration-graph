"""k50 — Exponentielle Ebene, Ecke abgeschnitten (BERICHT §29): gibt es eine geschlossene Kurve 𝒞 um 0 in D̄ ∩ f⁻¹(D̄), die L, L̄ im Abstand δ
umgeht, auf der der datengewichtete adjungierte Kopf exponentiell fällt?
  G_n(t) := (1/n)·ln|Λ_n(t)/Λ_n(L)| + ln(|f(t)|/ρ)      (f kubisches Modell, Λ_n adjungierter Kopf, ℓ_n = J̃_n^{−T}e_n)
Für jeden Winkel θ (Strahl von 0) wird min über r ∈ [r_min, ρ] von G_n(r e^{iθ}) mit |f(r e^{iθ})| ≤ ρ ausgegeben (Lage r*), dazu G_n auf Bögen um L im
Abstand δ (innen). Zusätzlich primal: (1/n)·ln max|δ_n/δ_n(L)| auf dem Rand von K_{r0}.
Aufruf: python k50_ecke_abschneiden.py 400
"""
import math
import os
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, acb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    n = int(sys.argv[1])
    t0 = time.time()
    bits = int(1.3 * n) + 400
    J = J_tilde("kubisch", n, bits=bits)
    ctx.prec = bits
    Jc = acb_mat(J)
    e = acb_mat(n + 1, 1); e[n, 0] = acb(1)
    c = Jc.solve(e, algorithm="approx")
    l = Jc.transpose().solve(e, algorithm="approx")
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    L = acb(arb(L_RE), arb(L_IM))
    Lam = acb_poly([l[k, 0] for k in range(n + 1)])      # in w = ρ/t
    Del = acb_poly([c[k, 0] for k in range(n + 1)])      # in z/ρ
    fk = acb_poly([acb(x) for x in kubik_koeff()])
    LamL = Lam(rho / L)
    DelL = Del(L / rho)
    rhof = float(rho.mid())

    def G(t):
        v = Lam(rho / t) / LamL
        a = abs(v)
        if not (a > 0):
            return None
        fv = abs(fk(t))
        return float((a.log() / n).mid()) + math.log(float(fv.mid()) / rhof), float(fv.mid()) / rhof

    print("n = %d (Lösen %.0f s)" % (n, time.time() - t0), flush=True)
    # Strahlen
    Lf = complex(float(L.real.mid()), float(L.imag.mid()))
    argL = math.degrees(math.atan2(Lf.imag, Lf.real))
    out = []
    for th in range(0, 360, 4):
        best = None
        for i in range(41):
            r = rhof * (0.40 + 0.60 * i / 40)
            t = acb(r * math.cos(math.radians(th)), r * math.sin(math.radians(th)))
            g = G(t)
            if g is None or g[1] > 1.0:
                continue
            if best is None or g[0] < best[0]:
                best = (g[0], r / rhof, g[1])
        out.append((th, best))
    for th, b in out:
        if b is None:
            print("  θ = %3d°: kein Punkt mit |f| ≤ ρ" % th)
        else:
            print("  θ = %3d°: min G = %+.4f bei r/ρ = %.3f (|f|/ρ = %.3f)" % (th, b[0], b[1], b[2]), flush=True)
    # Bögen um L im Abstand δ (Winkel von L aus)
    for dl in (0.003, 0.01, 0.03, 0.1):
        row = []
        for ph in range(130, 320, 15):
            t = acb(Lf.real + dl * math.cos(math.radians(ph)), Lf.imag + dl * math.sin(math.radians(ph)))
            g = G(t)
            row.append("%d°:%s" % (ph, "—" if g is None else ("%+.4f" % g[0]) + ("" if g[1] <= 1 else "(|f|>ρ)")))
        print("  δ = %.3f um L: " % dl + " ".join(row), flush=True)
    # primal auf dem Rand von K_{r0}
    for r0 in (0.1, 0.03):
        mx = -1e9
        for j in range(720):
            ph = 2 * math.pi * j / 720
            z = complex(rhof * math.cos(ph), rhof * math.sin(ph))
            if abs(z - Lf) < r0 or abs(z - Lf.conjugate()) < r0:
                continue
            v = abs(Del(acb(z.real, z.imag) / rho) / DelL)
            mx = max(mx, float((v.log() / n).mid()))
        print("  primal: max (1/n)ln|δ_n/δ_n(L)| auf ∂D ∖ U_%.2f = %+.5f" % (r0, mx), flush=True)
    print("fertig (%.0f s)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
