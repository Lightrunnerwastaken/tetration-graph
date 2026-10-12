"""Strang V — Antwort der Inversen auf KOHÄRENTE Störungen (Kneser, Tangente): η = J_N⁻¹T_N g mit glattem g (Newton-/Homotopie-Richtungen).

Hintergrund: v21 zeigt r_k ≈ e^{N g(k/N)} für Einzelgrade (g(1) ≈ 0.16) — keine diagonale Gewichtsnorm erfüllt zugleich (a) Lip polynomial und
(b) p_N − f_K klein. Die Familie in (H) besteht aber nur aus Fast-Lösungen (ARCHITEKTUR §5.1): Richtungen der Form J⁻¹(glatt).
Kandidat für die Norm: Graphennorm ‖η‖_J := ‖T_N J_∞ η‖_{ℓ¹_ρ} (Grösse des Residuums). Hier gemessen:
  Lip_J = ‖∂_ε (Ĵ^{f+εη})⁻¹‖ / ‖T_N g‖_{ℓ¹_ρ}  und  Lip_1 = dasselbe / ‖η‖_{ℓ¹_ρ}
für g ∈ {„kopf“: (x/ρ)^N (die (H-iii)-Richtung), „glatt“: Σ_{j≤20} r_j (x/ρ)^j/Σ|r_j|, „mitte“: Zufall in Graden 0.3N…0.7N, „eins“: 1,
„e<k>“: (x/ρ)^k (Ecken der ℓ¹-Kugel der Graphennorm), „echt“: η = p_N − T_N f_K (dann b = J̃ξ, Bedingung (b))}.
η = J̃_N⁻¹ b in der x-Koordinate (J̃_N = Ĵ mit k = ρw, ρ-skaliert, wie v14), m_P = m + ε(η'∘f + (m'/f')η), Linse S wie v12.
Aufruf: python v23_kohaerent.py N S art1,art2,...
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

ARGS = sys.argv[1:]
N = int(ARGS[0])
S = ARGS[1]
ARTEN = ARGS[2].split(",")
sys.argv = ["v20", str(N), S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12

sys.stdout.reconfigure(encoding="utf-8")
EPS = arb("1e-20")


def main():
    F, it, d = V12.theodorsen()
    bits = int(1.2 * N) + 400
    ctx.prec = bits
    rho = abs(V12.L_arb())
    a, m = V20.koeff_fK(N)
    fpol = arb_poly(a)
    mpol = arb_poly(m)
    fpr = fpol.derivative()
    mp = mpol.derivative()
    q = [arb(0)] * (N + 1)
    fpc = [fpr[i] if i <= fpr.degree() else arb(0) for i in range(N + 1)]
    mpc = [mp[i] if i <= mp.degree() else arb(0) for i in range(N + 1)]
    for i in range(N + 1):
        s_ = mpc[i]
        for j in range(i):
            s_ -= q[j] * fpc[i - j]
        q[i] = s_ / fpc[0]
    qpol = arb_poly(q)
    V12.modell_koeff = lambda n_: (a, m)
    Jx = V12.J_hat(N, np.zeros(N + 1, dtype=complex), bits)        # J̃_N (x-Koordinate, ρ-skaliert)
    ctx.prec = bits
    J0 = V12.J_hat(N, F, bits)
    ctx.prec = bits
    J0i = V20.inverse(J0)
    print(f"Kneser N = {N}, Linse S = {S}: ‖Ĵ⁻¹‖ = {np.linalg.norm(V20.zu_np(J0i), 2):.6f}", flush=True)
    rng = np.random.default_rng(7)
    for art in ARTEN:
        t0 = time.time()
        ctx.prec = bits
        bvec = [0.0] * (N + 1)                       # b = Koeffizienten von g in der Basis (x/ρ)^j
        if art == "kopf":
            bvec[N] = 1.0
        elif art == "eins":
            bvec[0] = 1.0
        elif art == "glatt":
            r = rng.standard_normal(21)
            for j in range(21):
                bvec[j] = r[j] / np.sum(np.abs(r))
        elif art == "mitte":
            lo, hi = int(0.3 * N), int(0.7 * N)
            r = rng.standard_normal(hi - lo)
            for j in range(lo, hi):
                bvec[j] = r[j - lo] / np.sum(np.abs(r))
        if art.startswith("e") and art[1:].isdigit():
            bvec[min(N, int(art[1:]))] = 1.0
        if art == "echt":                             # η = p_N − T_N f_K (echt), b = J̃ξ
            pn = [arb(x) for x in V20.PN[str(N)]]
            xi = arb_mat(N + 1, 1)
            for j in range(N + 1):
                xi[j, 0] = (pn[j] - a[j]) * rho ** j
            Bv = Jx * xi
            bvec = [float(Bv[j, 0].mid()) for j in range(N + 1)]
        else:
            B = arb_mat(N + 1, 1)
            for j in range(N + 1):
                B[j, 0] = arb(bvec[j])
            xi = Jx.solve(B, algorithm="approx")         # ξ_j = η_j ρ^j
        eta = [xi[j, 0] / rho ** j for j in range(N + 1)]
        n_eta = sum(abs(float(xi[j, 0].mid())) for j in range(N + 1))
        n_b = float(np.sum(np.abs(bvec)))
        top = sum(abs(float(xi[j, 0].mid())) for j in range(int(0.8 * N), N + 1))
        epol = arb_poly(eta)
        dep = epol.derivative()
        comp = arb_poly([0])                          # η'∘f mod x^{N+1} (Horner)
        for c in reversed(dep.coeffs()):
            comp = (comp * fpol).truncate(N + 1) + c
        dm = comp + (qpol * epol).truncate(N + 1)
        a1 = [a[j] + EPS * eta[j] for j in range(N + 1)]
        m1 = [m[j] + EPS * (dm[j] if j <= dm.degree() else 0) for j in range(N + 1)]
        V12.modell_koeff = lambda n_, a1=a1, m1=m1: (a1, m1)
        J1 = V12.J_hat(N, F, bits)
        ctx.prec = bits
        J1i = V20.inverse(J1)
        dn = np.linalg.norm(V20.zu_np(J1i - J0i), 2) / float(EPS.mid())
        print(f"N = {N:5d}  {art:6s}  ‖η‖_ℓ¹ρ = {n_eta:.3e} (Band ≥ 0.8N: {top:.2e})  ‖b‖ = {n_b:.2e}  Δ_inv/ε = {dn:.3e}  "
              f"Lip_J = {dn / n_b:.3e}  Lip_1 = {dn / n_eta:.3e}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
