"""O36 — Innenschicht-Profile der Eckmoden (Doppel-Skalenlimes an L): für die r grössten Singulärrichtungen von J_N⁻¹
(Eingang v_i, J_N v_i = u_i/s_i) speichern: Koeffizienten v_k (ℓ²_ρ), Bild J_∞v (Grade ≤ 3N), Kreisprofil |v(ρe^{iφ})| in
y = N(θ_L − φ). Ziel: Grenzgestalt in κ = k/N bzw. y prüfen (N = 200…600), Überlaufprofil j ↦ [J_∞v]_{N+j}.
Aufruf: python o36_innenschicht.py N [modell] [r]   → o36_<modell>_<N>.npz
"""
import sys

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
name = sys.argv[2] if len(sys.argv) > 2 else "kubisch"
r = int(sys.argv[3]) if len(sys.argv) > 3 else 3
ctx.prec = int(3.0 * N) + 400
NR = 3 * N if name != "kneser" else 3 * N
fk, mk, _ = K.modell(name, NR)
f = arb_poly(fk)
rho = (arb(K.L_RE) ** 2 + arb(K.L_IM) ** 2).sqrt()
pw = [arb_poly([1])]
for _ in range(N):
    pw.append((pw[-1] * f).truncate(NR + 1))
rp = [rho ** k for k in range(NR + 1)]
A = arb_mat(NR + 1, N + 1)
for j in range(N + 1):
    col = pw[j]
    for k in range(NR + 1):
        v = col[k] if k <= col.degree() else arb(0)
        if k >= j and k - j < len(mk):
            v = v + mk[k - j]
        A[k, j] = v * rp[k] / rp[j]
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
I = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    I[i, i] = arb(1)
Yi = JN.solve(I, algorithm="approx")
Y = np.array([[float(Yi[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
U, s, Vt = np.linalg.svd(Y)
print(f"{name} N={N}: s = " + ", ".join(f"{x:.4f}" for x in s[:6]), flush=True)
V = np.zeros((r, N + 1))
JV = np.zeros((r, NR + 1))
for i in range(r):
    yv = arb_mat(N + 1, 1)
    for j in range(N + 1):
        yv[j, 0] = arb(float(Vt[i, j]))
    va = JN.solve(yv, algorithm="approx")          # J_N va = Vt[i]  ⇒ va = s_i·U_i
    out = A * va
    nrm = np.sqrt(sum(float(va[j, 0].mid()) ** 2 for j in range(N + 1)))
    V[i] = [float(va[j, 0].mid()) / nrm for j in range(N + 1)]
    JV[i] = [float(out[k, 0].mid()) / nrm for k in range(NR + 1)]
th = np.angle(K.Lc)
ys = np.linspace(-300, 100, 801)
prof = np.array([np.abs(np.polynomial.polynomial.polyval(np.exp(1j * (th - ys / N)), V[i])) for i in range(r)])
np.savez(f"o36_{name}_{N}.npz", s=s[:10], V=V, JV=JV, ys=ys, prof=prof, N=N)
for i in range(r):
    k = np.arange(N + 1)
    w = V[i] ** 2
    print(f"Mode {i}: s = {s[i]:.4f};  Masse in k/N ≥ 0.9: {w[k >= 0.9 * N].sum():.3f}, ≥ 0.8: {w[k >= 0.8 * N].sum():.3f};  "
          f"Spitze |v| auf dem Kreis bei y = {ys[np.argmax(prof[i])]:+.1f}", flush=True)
