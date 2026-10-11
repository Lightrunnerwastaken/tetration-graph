"""X27 — Wo verschenkt Lemma X8 Faktoren? Die BESTE Konstante des Hauptterms
   A*(N) := max_{δ∈P_N} ‖δ'∘f‖_{L²(∂D)} / (N·‖J_∞δ‖_{ℓ²_ρ})      (J_∞δ bis Zeile 6N — mit 3N ist das Problem falsch gestellt: Zeilen 3N…5N können für das Extremal gross sein)
und die Teilkonstanten von Lemma X8:  (K/ρ)/min|f'|,  max|m/f'|/(ρc₀) (verbessert gegenüber max|m|/(ρc₀ min|f'|)),  max|m'/f'|.
Rechnung: A = T_{6N}J_∞|_{P_N} (ℓ²_ρ, arb), AᵀA = LLᵀ (Cholesky in arb), D_ij = δ_j'(f(x_i)) für δ_j = (x/ρ)^j an 8N Punkten x_i ∈ ∂D,
A*² = λ_max((L⁻¹Dᵀ)(L⁻¹Dᵀ)ᴴ)/(8N·N²).  Aufruf: python x27_atheta_opt.py N1 N2 …
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, acb_mat, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
mp.mp.dps = 40
rs = abs(mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337)))
rho = float(rs)
import os
PRECF = float(os.environ.get("PRECF", "3"))


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n):
    g = arb_poly([1 / a[0]])
    k = 1
    while k < n:
        k = min(2 * k, n)
        g = mullow(g, arb_poly([2]) - mullow(a.truncate(k), g, k), k)
    return g


# Randwerte (float, skalierte Reihen) für f, f', f'' auf |x| = ρ
K = len(A_str)
b0 = np.array([float(mp.mpf(a) * rs ** k) for k, a in enumerate(A_str)])
b1 = np.array([float(k * mp.mpf(A_str[k]) * rs ** k) for k in range(1, K)]) / rho
b2 = np.array([float(k * (k - 1) * mp.mpf(A_str[k]) * rs ** k) for k in range(2, K)]) / rho ** 2
for N in (int(a) for a in sys.argv[1:]):
    t0 = time.time()
    Ms = 8 * N
    phi = 2 * np.pi * (np.arange(Ms) + 0.5) / Ms
    w = np.exp(1j * phi)
    x = rho * w
    fx = np.polynomial.polynomial.polyval(w, b0)
    f1 = np.polynomial.polynomial.polyval(w, b1)
    f2 = np.polynomial.polynomial.polyval(w, b2)
    m = np.exp(x) / f1
    m1 = m * (1 - f2 / f1)
    if N == int(sys.argv[1]):
        print(f"Randgrössen: min|f'| = {np.abs(f1).min():.4f}, max|m| = {np.abs(m).max():.4f}, max|m/f'| = {np.abs(m / f1).max():.4f}, "
              f"max|m'/f'| = {np.abs(m1 / f1).max():.4f}", flush=True)
    # A = T_{3N} J_∞ |P_N in ℓ²_ρ (arb)
    Mr = int(sys.argv[0] and 6) * N
    ctx.prec = int(PRECF * N) + 600
    F = arb_poly([arb(s) for s in A_str[: Mr + 1]])
    ex = arb_poly([1 / arb.fac_ui(k) for k in range(Mr + 1)])
    mm = mullow(ex, inv_series(F.derivative(), Mr + 1), Mr + 1)
    rk = [arb(1)]
    for k in range(1, Mr + 1):
        rk.append(rk[-1] * arb(mp.nstr(rs, 50)))
    A = arb_mat(Mr + 1, N + 1)
    pw = arb_poly([1])
    for j in range(N + 1):
        for k in range(Mr + 1):
            A[k, j] = (pw[k] + (mm[k - j] if k >= j else 0)) * rk[k] / rk[j]
        pw = mullow(pw, F, Mr + 1)
    G = A.transpose() * A
    # Cholesky in arb über mpmath-freie einfache Implementierung (G symmetrisch positiv definit)
    n = N + 1
    Lc = arb_mat(n, n)
    for j in range(n):
        s = G[j, j]
        for k in range(j):
            s -= Lc[j, k] * Lc[j, k]
        Lc[j, j] = arb(s.mid()).sqrt()
        Lc[j, j] = arb(Lc[j, j].mid())
        for i in range(j + 1, n):
            s = G[i, j]
            for k in range(j):
                s -= Lc[i, k] * Lc[j, k]
            Lc[i, j] = arb((s / Lc[j, j]).mid())
    # D^T: Zeile j = Werte von δ_j'(f(x_i)) = j f(x_i)^{j−1}/ρ^j  (komplex) → Real- und Imaginärteil getrennt
    Dt_re = arb_mat(n, Ms)
    Dt_im = arb_mat(n, Ms)
    # f(x_i) in acb aus der (abgeschnittenen) Reihe — alle Potenzen exakt aus DIESEM Wert (konsistent; float verliert 1.62^N)
    rho_a = arb(mp.nstr(rs, 60))
    koeff = [arb(s) for s in A_str[:2500]]
    fxa = []
    for i in range(Ms):
        xi = acb(rho_a * arb(phi[i]).cos(), rho_a * arb(phi[i]).sin())
        v = acb(0)
        for c in reversed(koeff):
            v = v * xi + c
        fxa.append(v)
    pw_a = [acb(1)] * Ms
    rpow = arb(1)
    for j in range(n):
        if j >= 1:
            rpow = rpow * rho_a
            for i in range(Ms):
                val = j * pw_a[i] / rpow
                Dt_re[j, i] = arb(val.real.mid())
                Dt_im[j, i] = arb(val.imag.mid())
            pw_a = [pw_a[i] * fxa[i] for i in range(Ms)]
    def vorwaerts(B):
        X = arb_mat(n, B.ncols())
        for i in range(n):
            for c in range(B.ncols()):
                s = B[i, c]
                for k in range(i):
                    s -= Lc[i, k] * X[k, c]
                X[i, c] = arb((s / Lc[i, i]).mid())
        return X
    Xre = vorwaerts(Dt_re)      # L⁻¹ Dᵀ, Vorwärtseinsetzen auf Mittelpunkten (Schätzung, nicht streng)
    Xim = vorwaerts(Dt_im)
    Xr = np.array([[float(Xre[i, j].mid()) for j in range(Ms)] for i in range(n)])
    Xi = np.array([[float(Xim[i, j].mid()) for j in range(Ms)] for i in range(n)])
    Xc = Xr + 1j * Xi
    smax = np.linalg.svd(Xc, compute_uv=False)[0]
    Astern = smax / np.sqrt(Ms) / N
    print(f"N = {N}: A*(N) = max ‖δ'∘f‖/(N‖J_∞δ‖) = {Astern:.4f}  ({time.time() - t0:.0f}s)", flush=True)
