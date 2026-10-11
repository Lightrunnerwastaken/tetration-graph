"""D30 — Rest-Schranke (a) über die Faktorisierung J_N = M_N (I + K_N) (kubisches Modell).

M_N = P_N M_m P_N (untere Dreiecks-Toeplitz von m = f'∘f, exakt invertierbar, M_N⁻¹ = P_N M_{1/m} P_N),
K_N = M_N⁻¹ C_N = P_N K P_N, Kδ = (δ∘f)/(f'∘f).  Dann J_N⁻¹ = (I + K_N)⁻¹ M_N⁻¹ und
    σ_{r+1}(J_N⁻¹) ≤ σ_{r+1}((I + K_N)⁻¹) · ‖M_N⁻¹‖ ≤ σ_{r+1}((I + K_N)⁻¹) · max_{|x|=|L|} |1/m(x)|.
Gemessen: max_{∂D}|1/m| gegen max_{Γin}|1/m| = ‖σ‖_∞; Singulärwerte von Y := (I + K_N)⁻¹ = J_N⁻¹M_N (ℓ²_{|L|});
Vergleich mit dem Toeplitz-Abschnitt T_N(χ), χ = 1_{Γin} (Parametrix-Kandidat): Singulärwerte von Y − T_N(χ).
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
L = complex(0.31813150520476413, 1.3372357014306895)
R = abs(L)
mu = np.sqrt(L)
rL = (mu - 1) / (L - L.conjugate())
beta = rL.imag / L.imag
alpha = rL.real - beta * L.real
coef = [R ** 2 * alpha, 1 + R ** 2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]
f = np.polynomial.Polynomial(coef)
fp = f.deriv()
m = fp(f)
Mq = 1 << 15
phi = 2 * np.pi * np.arange(Mq) / Mq
x = R * np.exp(1j * phi)
inn = np.abs(f(x)) < R
im = np.abs(1 / m(x))
print(f"max_(∂D)|1/m| = {im.max():.5f} bei {np.degrees(phi[im.argmax()]):.2f}°;  max_(Γin)|1/m| = {im[inn].max():.5f};  "
      f"max_(Γout)|1/m| = {im[~inn].max():.5f}", flush=True)
chi_c = np.fft.fft(inn.astype(float)) / Mq          # Fourier-Koeff. von χ: χ_d = chi_c[d % Mq]
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(0.9 * N) + 300
    a = [arb(float(c)) for c in coef] + [arb(0)] * (n - 4)
    _, J = K.FJ(a, N)
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    sc = R ** k
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * sc[:, None] / sc[None, :]
    mc = np.zeros(n)
    mco = m.coef
    mc[:min(n, len(mco))] = mco[:min(n, len(mco))]
    Mt = np.zeros((n, n))
    for d in range(min(n, len(mco))):
        Mt += np.diag(np.full(n - d, mc[d] * R ** d), -d)        # umskaliert: m_d R^d auf der d-ten Unterdiagonale
    Y = G @ Mt
    sY = np.linalg.svd(Y, compute_uv=False)
    sG = np.linalg.svd(G, compute_uv=False)
    T = np.array([[chi_c[(i - j) % Mq].real for j in range(n)] for i in range(n)])
    sT = np.linalg.svd(T, compute_uv=False)
    sD = np.linalg.svd(Y - T, compute_uv=False)
    print(f"N={N} ({time.time()-t0:.0f}s): σ(J⁻¹) = {sG[0]:.3f}, {sG[1]:.3f}, {sG[2]:.4f}, {sG[3]:.4f} | σ((I+K)⁻¹) = {sY[0]:.3f}, "
          f"{sY[1]:.3f}, {sY[2]:.4f}, {sY[3]:.4f}, {sY[4]:.4f} | ‖T_N(χ)‖ = {sT[0]:.4f} | σ(Y − T_N(χ)) = {sD[0]:.3f}, {sD[1]:.3f}, "
          f"{sD[2]:.4f}, {sD[3]:.4f}, {sD[5]:.4f}, {sD[10]:.4f}", flush=True)
