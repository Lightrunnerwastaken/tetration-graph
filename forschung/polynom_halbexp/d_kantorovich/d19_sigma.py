"""D19 — Wachstum der Eckmode bis N = 1000, billig: nur die letzten B Spalten von J_N⁻¹ (J bei Π_N f_K, Strang A),
σ1, σ2 von |L|-umskaliertem G[:, N−B:], Lage der Ausgabe-Spitze, und die Antwort auf die echte Eingabe
g = J∞ f_K = e^x + (f_K'∘f_K) f_K (über die Identität J_N Π_N f_K = P_N g: Komponente von Π_N f_K entlang u1).
Aufruf: python d19_sigma.py N1,N2,... [B] [prec_faktor]
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
B = int(sys.argv[2]) if len(sys.argv) > 2 else 60
PF = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
PI = json.load(open(K.GEMEINSAM + "/lin_PiN_fK_bis_1000.json"))
L = 1.3745570107437075
datei = K.HIER + "/d19_sigma.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(PF * N) + 300
    a = [arb(s) for s in PI[str(N)]["b"][:n]]
    _, J = K.FJ(a, N)
    E = arb_mat(n, B)
    for c in range(B):
        E[N - B + 1 + c, c] = 1
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * (L ** k)[:, None] / (L ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    M = 8192
    Fu = np.abs(np.fft.ifft(np.concatenate([U[:, 0], np.zeros(M - n)])) * M)
    # Komponente der echten Lösung Π_N f_K entlang u1 (umskaliert): zeigt, ob die Eingabe die Eckmode anregt
    x = np.array([float(arb(s).mid()) for s in PI[str(N)]["b"][:n]]) * L ** k
    z = {"s1": float(S[0]), "s2": float(S[1]), "s3": float(S[2]), "u1_phi": 360.0 * int(np.argmax(Fu)) / M,
         "Pi_entlang_u1": float(abs(U[:, 0] @ x)), "Pi_entlang_u2": float(abs(U[:, 1] @ x)),
         "v1_letzte": [float(abs(Vt[0, B - 1 - d])) for d in (0, 2, 5, 10, 20, 30, 40, 50) if d < B]}
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    print(f"N={N} ({time.time()-t0:.0f}s, prec {ctx.prec}) σ1={S[0]:.4f} σ2={S[1]:.4f} σ3={S[2]:.4f} | u1 bei {z['u1_phi']:.2f}° | "
          f"Π_N f_K·u1 = {z['Pi_entlang_u1']:.2e}, ·u2 = {z['Pi_entlang_u2']:.2e} | |v1| ab Spitze: {[f'{t:.0e}' for t in z['v1_letzte']]}",
          flush=True)
