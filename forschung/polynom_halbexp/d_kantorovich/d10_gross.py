"""D10 — Stabilitätsnormen und Eckspitze für grosse N (bis 250): |J_N^{-1}| in ℓ1/ℓ2 mit w_k = |L|^k(1+k)^s,
Eckspitze H_N (sup auf |x| = |L| der Eckspalten), Eckblock-Normen. Fortsetzung von d07/d09.
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
p = K.loesungen(max(NS), bits=1200)
L = 1.3745570107437075
datei = K.HIER + "/d10_gross.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
t0 = time.time()
for N in NS:
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    _, J = K.FJ(a, N)
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    z = {}
    for s in (0, 1, 2, 3, 4):
        w = (1.0 + k) ** s
        A = G * w[:, None] / w[None, :]
        col = np.abs(A).sum(axis=0)
        z[f"l1_{s}"] = float(col.max())
        z[f"l1j_{s}"] = int(col.argmax())
        z[f"l2_{s}"] = float(np.linalg.norm(A, 2))
        z[f"l2ecke_{s}"] = float(np.linalg.norm(A[:, N - 30:], 2))
        z[f"l2rest_{s}"] = float(np.linalg.norm(A[:, :N - 30], 2))
    M = 8192
    best = (0.0, None, None)
    for j in range(N - 40, n):
        F = np.abs(np.fft.ifft(np.concatenate([G[:, j], np.zeros(M - n)])) * M)
        i = int(np.argmax(F))
        if F[i] > best[0]:
            best = (float(F[i]), j, 360.0 * i / M)
    z["H"], z["H_j"], z["H_phi"] = best
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    print(f"N={N} ({time.time()-t0:.0f}s) H={z['H']:.3f}@{z['H_phi']:.2f}° (j={z['H_j']}) | "
          + " ".join(f"s={s}: l1={z[f'l1_{s}']:.3f} l2={z[f'l2_{s}']:.3f} [Ecke {z[f'l2ecke_{s}']:.3f} Rest {z[f'l2rest_{s}']:.3f}]"
                     for s in (0, 1, 2, 3, 4)), flush=True)
