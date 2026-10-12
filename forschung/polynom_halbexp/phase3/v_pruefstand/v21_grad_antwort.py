"""Strang V — Antwort der Inversen auf Störungen EINES Grades (Kneser, Tangente): r_k(N) = ‖∂_ε (Ĵ_N^{f_K + ε(x/ρ)^k})⁻¹‖ bei ε = 0.

Zweck (Auftrag Leitung 22:52): Gewichtete Norm ‖η‖_ω = Σ_k |η_k| ρ^k ω_k für Störungen δP suchen. Erster Ordnung ist die Lipschitz-Konstante
von P ↦ (Ĵ^P)⁻¹ in dieser Norm genau sup_k r_k(N)/ω_k (ℓ¹-Kugel: Ecken = Einzelgrade). Also: r_k messen, dann ω wählen.
P = f_K + εη, η = (x/ρ)^k; m_P = P'∘P = m + ε(η'∘f + (f''∘f)η) + O(ε²), f''∘f = m'/f' (aus m = f'∘f), η'∘f = k f^{k−1}/ρ^k.
Ĵ wie v12/v20 (konforme Abbildung, Linse S, Ks Profilform), f_K aus gemeinsam/kneser_halbexp_koeff_4000.json.
r_k = ‖(Ĵ^{P})⁻¹ − Ĵ⁻¹‖/ε mit ε = 1e-20 (Differenzenquotient; Fehler O(ε)).
Aufruf: python v21_grad_antwort.py N S [k1,k2,...|auto]
"""
import sys
import time

import numpy as np
from flint import arb, arb_poly, ctx

ARGS = sys.argv[1:]
N = int(ARGS[0])
S = ARGS[1]
KS = ARGS[2] if len(ARGS) > 2 else "auto"
sys.argv = ["v20", str(N), S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12

sys.stdout.reconfigure(encoding="utf-8")
EPS = arb("1e-20")


def potenz(P, e, n):
    r = arb_poly([1])
    b = P
    while e > 0:
        if e & 1:
            r = (r * b).truncate(n + 1)
        e >>= 1
        if e:
            b = (b * b).truncate(n + 1)
    return r


def main():
    if KS == "auto":
        ks = sorted(set([1, 2, 5, 10, 20, 50] + [int(round(x * N)) for x in (0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98)] + [N]))
    else:
        ks = [int(x) for x in KS.split(",")]
    F, it, d = V12.theodorsen()
    bits = int(1.2 * N) + 400
    ctx.prec = bits
    rho = abs(V12.L_arb())
    a, m = V20.koeff_fK(N)
    fpol = arb_poly(a)
    mpol = arb_poly(m)
    fpr = fpol.derivative()
    # f''∘f = m'/f' als Reihe (Division durch f', f'(0) ≠ 0)
    mp = mpol.derivative()
    q = [arb(0)] * (N + 1)
    fpc = [fpr[i] if i <= fpr.degree() else arb(0) for i in range(N + 1)]
    mpc = [mp[i] if i <= mp.degree() else arb(0) for i in range(N + 1)]
    for i in range(N + 1):
        s_ = mpc[i]
        for j in range(i):
            s_ -= q[j] * fpc[i - j]
        q[i] = s_ / fpc[0]
    V12.modell_koeff = lambda n_: (a, m)
    J0 = V12.J_hat(N, F, bits)
    ctx.prec = bits
    J0i = V20.inverse(J0)
    X0 = V20.zu_np(J0i)
    print(f"Kneser N = {N}, Linse S = {S}: ‖Ĵ⁻¹‖ = {np.linalg.norm(X0, 2):.6f}", flush=True)
    for k in ks:
        t0 = time.time()
        ctx.prec = bits
        rk = rho ** k
        a1 = list(a)
        a1[k] = a1[k] + EPS / rk
        fk1 = potenz(fpol, k - 1, N) if k >= 1 else arb_poly([0])
        dm = (fk1 * (arb(k) / rk)) if k >= 1 else arb_poly([0])
        qx = (arb_poly(q) * arb_poly([0] * k + [1 / rk])).truncate(N + 1)
        dm = dm + qx
        m1 = [m[i] + EPS * (dm[i] if i <= dm.degree() else 0) for i in range(N + 1)]
        V12.modell_koeff = lambda n_, a1=a1, m1=m1: (a1, m1)
        J1 = V12.J_hat(N, F, bits)
        ctx.prec = bits
        J1i = V20.inverse(J1)
        D = V20.zu_np(J1i - J0i)
        r = np.linalg.norm(D, 2) / float(EPS.mid())
        print(f"N = {N:5d}  k = {k:5d}  ξ = {k / N:.3f}  r_k = {r:.4e}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
