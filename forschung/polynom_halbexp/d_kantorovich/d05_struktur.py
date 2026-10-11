"""D05 — Struktur von J_N und J_N^{-1} in den Gewichten w_k = |L|^k (1+k)^s: Wärmebilder, Spaltennormen, Ecke.

Zerlegung (exakt, Lemma D3): J_N = U_N S_N + M_N, U = Komposition mit u = p − a_0 (untere Dreiecksmatrix),
S = Taylor-Verschiebung nach a_0 (obere Dreiecksmatrix), M = Multiplikation mit p'∘p (untere Dreiecks-Toeplitz).
"""
import sys

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
p = K.loesungen(N, bits=1200)
ctx.prec = 500 + 3 * N
absL = 1.3745570107437075
a = [arb(x.mid()) for x in p[N]]
P = K.potenzen(a, N)
_, J = K.FJ(a, N, P)
n = N + 1
Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
Jf = np.array([[float(J[i, j].mid()) for j in range(n)] for i in range(n)])
Jif = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)])
# Komposition C_kj = [x^k] p^j und Multiplikation M
Cf = np.array([[float(P[j][i].mid()) for j in range(n)] for i in range(n)])
Mf = Jf - Cf
k = np.arange(n)
fig, ax = plt.subplots(2, 3, figsize=(15, 9))
for row, s in enumerate((0, 4)):
    w = absL ** k * (1.0 + k) ** s
    for col, (M, name) in enumerate(((Cf, "C (Komposition)"), (Mf, "M (Multiplikation)"), (Jif, "J^{-1}"))):
        R = np.abs(M) * w[:, None] / w[None, :]
        im = ax[row, col].imshow(np.log10(R + 1e-300), vmin=-12, vmax=8, cmap="viridis")
        ax[row, col].set_title(f"log10 |{name}| gewichtet, s={s}")
        plt.colorbar(im, ax=ax[row, col])
    colnorm = (np.abs(Jif) * w[:, None]).sum(axis=0) / w
    print(f"s={s}: Spaltennormen von J^-1 (gewichtet), j=0,1,2,5,10,N/4,N/2,3N/4,N-5,N:")
    print("   ", [round(float(colnorm[j]), 3) for j in (0, 1, 2, 5, 10, N // 4, N // 2, 3 * N // 4, N - 5, N)])
    rownorm = (np.abs(Jif) * w[:, None]).sum(axis=1)
    # wo sitzt die Masse von J^{-1} e_0?
    v0 = np.abs(Jif[:, 0]) * w
    print("    Masse von J^-1 e_0 auf k: Anteile k<10, 10..N/2, >N/2:", [round(float(x), 3) for x in (v0[:10].sum(), v0[10:N // 2].sum(), v0[N // 2:].sum())])
    vN = np.abs(Jif[:, N]) * w / w[N]
    print("    Masse von J^-1 e_N/w_N: k<10, 10..N/2, >N/2:", [float(f"{x:.3e}") for x in (vN[:10].sum(), vN[10:N // 2].sum(), vN[N // 2:].sum())])
plt.tight_layout()
plt.savefig(K.HIER + f"/d05_struktur_N{N}.png", dpi=80)
# Diagonaldominanz der Komposition in der Ecke
w = absL ** k
Cd = np.abs(np.diag(Cf)) ; Md = np.abs(np.diag(Mf))
print("Diagonale C_jj (j=10,N/2,N):", [f"{Cd[j]:.3e}" for j in (10, N // 2, N)], " M_jj:", f"{Md[N]:.4f}")
np.savez_compressed(K.HIER + f"/d05_N{N}.npz", J=Jf, Ji=Jif, C=Cf)
