"""W20 — Eckanteil des kopffreien Pol-Antwort T(η) := ⟨y, B_N[k_η]_{>N}⟩ (= ⟨u₁, y_η⟩ bei Cs y_η = J⁻¹B_N[k_η]_{>N}),
y = J^{-T}u₁ hochgenau, auf ∂W̄: A-Bogen und D-Bogen (Γin), kubisches Modell, unskaliert.
B_N[k_η]_{>N} = η^{−N−1}[f^{N+1}/(η − f)]_{≤N}.
Aufruf: python w20_eckanteil_rand.py N1,N2,...
"""
import sys
import time

import numpy as np
from flint import acb, acb_series, arb, arb_mat, arb_poly, acb_poly, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho = W.RHO
L = W.L
a0 = coef[0]
R0 = abs(L - a0)
psiL = np.angle(L - a0)
argL = np.angle(L)
pts = []
for dist in (0.0, 0.02, 0.05, 0.1, 0.19, 0.38, 0.6, 1.0):            # A-Bogen, Abstand von L
    dpsi = 2 * np.arcsin(dist / (2 * R0))
    pts.append((f"A{dist}", a0 + R0 * np.exp(1j * (psiL - dpsi))))
for tau in (2, 8, 24, 103.38):                                         # D-Bogen
    pts.append((f"D{tau}", rho * np.exp(1j * (argL + np.deg2rad(tau)))))
for N in [int(x) for x in sys.argv[1].split(",")]:
    n = N + 1
    ctx.prec = int(1.5 * N) + 400
    ctx.cap = N + 10
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    u1 = d["U"][:, 0]
    t0 = time.time()
    a = W.kubik_a(N, coef)
    J = W.J_mat(a, N)
    rr = arb(rho)
    Dk = [rr ** k for k in range(n)]
    Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)]), algorithm="approx")
    # unskaliert: ℓ(g) = Σ_j (y_j / ρ^j) g_j  mit y_j = Z_j/ρ^j·... :  J̃ᵀy = u, J̃ = DJD⁻¹ ⇒ Z = D y; Paarung ⟨y, D g⟩ = Σ Z_j/ρ^j ρ^j g_j/ρ^j…
    # sorgfältig: J̃ᵀ y = u ⇔ Jᵀ(Dy) = Du. Eckanteil der unskalierten Antwort: ⟨u, D J⁻¹ g⟩ = ⟨J^{-T} D u, g⟩ = ⟨Z, g⟩.
    fpol = acb_poly([acb(float(c)) for c in coef])
    fp = fpol.derivative()
    mser = acb_poly([0]); Pw = acb_poly([1])
    for c in fp.coeffs():
        mser += c * Pw; Pw = (Pw * fpol).truncate(n)
    mser = mser.truncate(n)
    uD = [arb(float(u1[k])) * Dk[k] for k in range(n)]
    out = []
    for name, eta in pts:
        et = acb(eta.real, eta.imag)
        s1 = list(acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs())
        kz = [et ** (-(j + 1)) for j in range(n)]
        s2 = list((mser * acb_poly(kz)).truncate(n).coeffs())
        s1 += [acb(0)] * (n - len(s1)); s2 += [acb(0)] * (n - len(s2))
        s = acb(0)
        for j in range(n):
            s += Z[j, 0] * (s1[j] + s2[j]) - uD[j] * kz[j]
        m = abs(s).mid()
        out.append((name, float(m.log()) / N if m > 0 else -1e9))
    print(f"N={N} ({time.time()-t0:.0f}s): " + " ".join(f"{nm}:{v:+.4f}" for nm, v in out), flush=True)
