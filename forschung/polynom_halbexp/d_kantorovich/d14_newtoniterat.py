"""D14 — Stabilität für grosse N beim ersten Newton-Iterat q1 = q_N − J(q_N)^{-1}F(q_N) (liegt bis auf ~1e-12 an p_N,
§4.7), ohne p_N selbst zu berechnen. Ausgabe: |J(q1)^{-1}| in ℓ1/ℓ2 (w_k = |L|^k (1+k)^s), Eckblock, Eckspitze H_N.
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
kn = json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_1000.json"))["a"]
L = 1.3745570107437075
datei = K.HIER + "/d14_newtoniterat.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
t0 = time.time()
for N in NS:
    ctx.prec = 500 + 3 * N
    q = [arb(str(x)) for x in kn[:N + 1]]
    F, J = K.FJ(q, N)
    d = K.loese(J, F)
    q1 = [arb((q[i] - d[i]).mid()) for i in range(N + 1)]
    F1, J1 = K.FJ(q1, N)
    d1 = K.loese(J1, F1)
    n = N + 1
    k = np.arange(n)
    schritt2 = float(sum(abs(float(d1[i].mid())) * L ** i * (1 + i) ** 4 for i in range(n)))
    Ji = J1.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    z = {"zweiter_newtonschritt_l1_4": schritt2}
    for s in (0, 2, 3, 4):
        w = (1.0 + k) ** s
        A = G * w[:, None] / w[None, :]
        z[f"l1_{s}"] = float(np.abs(A).sum(axis=0).max())
        z[f"l2_{s}"] = float(np.linalg.norm(A, 2))
        z[f"l2ecke_{s}"] = float(np.linalg.norm(A[:, N - 30:], 2))
        z[f"l2rest_{s}"] = float(np.linalg.norm(A[:, :N - 30], 2))
    M = 8192
    best = (0.0, None, None)
    for j in range(N - 40, n):
        Fj = np.abs(np.fft.ifft(np.concatenate([G[:, j], np.zeros(M - n)])) * M)
        i = int(np.argmax(Fj))
        if Fj[i] > best[0]:
            best = (float(Fj[i]), j, 360.0 * i / M)
    z["H"], z["H_j"], z["H_phi"] = best
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    print(f"N={N} ({time.time()-t0:.0f}s) 2. Newton-Schritt {schritt2:.1e} | H={z['H']:.3f}@{z['H_phi']:.2f}° | "
          + " ".join(f"s={s}: l1={z[f'l1_{s}']:.3f} l2={z[f'l2_{s}']:.3f} [Ecke {z[f'l2ecke_{s}']:.3f} Rest {z[f'l2rest_{s}']:.3f}]"
                     for s in (0, 2, 3, 4)), flush=True)
