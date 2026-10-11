"""b30: Spitzen-Anteil der Pol-Antwort für Pole im Abstand u = v·N^{−0.883} von L (Winkelhalbierende).
Spitzenbeitrag := σ_1·|⟨v_1, w̃_ζ⟩| (gewichtet ℓ², w_k = |L|^k);  normale Antwort := |Ã w̃_ζ| ohne Spitze.
Wächst σ_1·Überlapp für feste v mit N? Dann würde (H_R) für Pole extrem nahe an L irgendwann von der Spitze dominiert."""
import sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
aL = abs(L); a0 = mp.mpf(K[0])
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
for N in (160, 200, 250, 300, 360):
    pw, g, Ji = baue(N, J=N)
    w = np.array([float(aL ** k) for k in range(N + 1)])
    Jf = np.array([[float(Ji[i, j].mid().str(30, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    A = (w[:, None] * Jf) / w[None, :]
    U, sv, Vt = np.linalg.svd(A)
    zeile = []
    for v in (0.3, 1.0, 4.0):
        u = mp.mpf(v) * mp.mpf(N) ** (-0.883)
        zeta = L + u * mp.expj(psi)
        rt = np.array([complex(x) for x in rhs_pol(N, pw, g, zeta)]) * w
        ov = abs(np.dot(Vt[0], rt))
        rest = rt - np.dot(Vt[0], rt) * Vt[0]
        normal = np.linalg.norm(A @ rest)
        zeile.append(f"v={v}: Spitze {sv[0]*ov:.3e} / normal {normal:.3e}")
    print(f"N={N} (σ1={sv[0]:.3f}, σ2={sv[1]:.3f}): " + " | ".join(zeile), flush=True)
