"""M06 — Pol-Antwort an der Ecke im normierten radialen Quadrat-Modell: A_ζ(z) := −(ζ − z)·E_N(ζ, z), E_N = Π_N k_ζ − k_ζ.
A_ζ ist ein Polynom vom Grad N+1 mit A_ζ(ζ) = 1 und T_N[A∘f·(ζ − z) + m·A·(ζ − f)] = 0.
Ausgabe: Koeffizientenprofil von A, |A| auf Kreisen, |E| auf dem Einheitskreis.
Aufruf: python m06_pol_ecke.py N eps theta_grad "zeta-Liste (Python-Ausdrücke, L = e^{iθ})"
"""
import sys
import time

import numpy as np
from flint import acb, acb_mat, arb, ctx

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
eps = float(sys.argv[2])
thg = float(sys.argv[3])
th = np.deg2rad(thg)
L = np.exp(1j * th)
zetas = [complex(eval(s, {"np": np, "L": L})) for s in sys.argv[4].split(";")]
coef = [eps, 1 - 2 * eps * np.cos(th), eps]
ctx.prec = int(1.2 * N) + 300
t0 = time.time()
J = K.jmat(coef, N, reell=True)
# komplexe Kopie
Jc = acb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        Jc[i, j] = acb(J[i, j])
print(f"J aufgebaut ({time.time()-t0:.0f}s)", flush=True)
# Taylor-Koeffizienten von J_∞ k_ζ bis Grad N: k_ζ∘f + m k_ζ
from flint import acb_poly
fpoly = acb_poly([acb(c) for c in coef])
c0 = 1 - 2 * eps * np.cos(th)
m = acb_poly([acb(c0), acb(0), acb(0)]) + acb(2 * eps) * fpoly    # m = c0 + 2ε f
for zeta in zetas:
    z = acb(zeta.real, zeta.imag)
    # k_ζ(x) = Σ x^j ζ^{-j-1}; k_ζ∘f = 1/(ζ − f): Reihe über Newton/Division
    den = acb_poly([z]) - fpoly
    kf = (acb_poly([1]) / den) if False else None
    # Potenzreihen-Inverse von (ζ − f) bis Grad N
    d = [den[i] for i in range(3)]
    inv = [acb(0)] * (N + 1)
    inv[0] = 1 / d[0]
    for n in range(1, N + 1):
        s = acb(0)
        for i in range(1, 3):
            if n - i >= 0:
                s += d[i] * inv[n - i]
        inv[n] = -s / d[0]
    kz = [z ** (-(j + 1)) for j in range(N + 1)]
    mk = [acb(0)] * (N + 1)
    for j in range(N + 1):
        for i in range(3):
            if j - i >= 0:
                mk[j] += m[i] * kz[j - i]
    rhs = acb_mat(N + 1, 1, [inv[j] + mk[j] for j in range(N + 1)])
    X = Jc.solve(rhs, algorithm="approx")
    p = np.array([complex(X[j, 0].mid()) for j in range(N + 1)])     # Π_N k_ζ
    # A = 1 − (ζ − x)·p(x)  (Grad N+1)
    A = np.zeros(N + 2, complex)
    A[0] = 1
    A[: N + 1] -= zeta * p
    A[1:] += p
    print(f"\nζ = {zeta:.6f} (|ζ| = {abs(zeta):.6f}) ({time.time()-t0:.0f}s)")
    a = np.abs(A)
    print("  |A_k|: " + ", ".join(f"k={k}:{a[k]:.2e}" for k in list(range(0, N + 2, max(1, N // 16))) + [N - 5, N, N + 1]))
    phis = np.deg2rad(np.array([0, 30, 60, 70, 74, 76, thg, 78, 80, 90, 120, 150, 180]))
    for r in (0.5, 0.8, 0.95, 1.0):
        vals = [abs(np.polyval(A[::-1], r * np.exp(1j * ph))) for ph in phis]
        print(f"  |A| auf r={r}: " + " ".join(f"{v:.2e}" for v in vals))
    # E auf dem Einheitskreis
    vals = [abs(np.polyval(p[::-1], np.exp(1j * ph)) - 1 / (zeta - np.exp(1j * ph))) for ph in phis]
    print("  |E| auf r=1:  " + " ".join(f"{v:.2e}" for v in vals))
    print("  Winkel:       " + " ".join(f"{np.rad2deg(ph):8.1f}" for ph in phis))
