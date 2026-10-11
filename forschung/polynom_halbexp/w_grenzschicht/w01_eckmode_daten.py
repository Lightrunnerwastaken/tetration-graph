"""W01 — Kubisches Modell: volle umskalierte Inverse, SVD, Eckmoden u_i, v_i speichern (N moderat).
Aufruf: python w01_eckmode_daten.py 200,300,400
Ausgabe: w01_N{N}.npz mit S (Singulärwerte), U (erste 4 linke = Ausgabe), Vt (erste 4 rechte = Eingabe), coef.
"""
import sys
import time

import numpy as np
from flint import ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
coef, Lf = W.kubik_koeff()
print("coef", coef, flush=True)
t0 = time.time()
for N in NS:
    ctx.prec = int(0.9 * N) + 300
    a = W.kubik_a(N, coef)
    J = W.J_mat(a, N)
    G = W.inverse_umskaliert(J, N)
    U, S, Vt = np.linalg.svd(G)
    np.savez(W.HIER + f"/w01_N{N}.npz", S=S, U=U[:, :6], Vt=Vt[:6], coef=np.array(coef), G=G)
    print(f"N={N} ({time.time()-t0:.0f}s): σ = {S[0]:.5f} {S[1]:.5f} {S[2]:.5f} {S[3]:.5f}", flush=True)
