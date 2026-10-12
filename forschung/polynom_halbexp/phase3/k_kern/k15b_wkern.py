"""k15b — Volterra-Kern in der w-Koordinate: H_n = ‖c_n‖ e_nᵀ/M_nn − Σ_j K̂_nj H_j mit M = C_k L (LU-Faktor von Ĵ), K̂_nj = ‖c_n‖M_nj/(‖c_j‖M_nn).
Frage: ist der w-Kern beschränkt (dann taugt die Rekursion für eine Kontraktion) oder ebenso exponentiell wie in x (k15)?"""
import sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import k04_konform as K4
from flint import ctx
sys.stdout.reconfigure(encoding="utf-8")
d = np.load("k15_volterra_kern_kubisch_300.npz")
L = d["L"]; h = d["h"]; q = d["q"]; Nr = L.shape[0] - 1
cn = q / h
ctx.prec = 800
F, _ = K4.theodorsen()
C = K4.float_mat(K4.C_matrix(F, Nr))
M = C @ L
for n in [100, 150, 200, 250, Nr]:
    row = np.array([cn[n] * M[n, j] / (cn[j] * M[n, n]) for j in range(n)])
    a = np.abs(row)
    print("n = %d: Σ_j|K̂_nj| = %.4g, Maximum bei j = %d (t = %.3f), |K̂_{n,n−1}| = %.3e, |K̂_{n,n−5}| = %.3e, Σ über j < 0.7n: %.3e"
          % (n, a.sum(), np.argmax(a), np.argmax(a) / n, a[n - 1], a[n - 5], a[: int(0.7 * n)].sum()))
