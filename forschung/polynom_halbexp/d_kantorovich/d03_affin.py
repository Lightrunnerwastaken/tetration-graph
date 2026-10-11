"""D03 — affin-invariantes Kantorovich (Deuflhard) in gewichteten ℓ1-Normen: Rohdaten.

Startpunkt x0: 'pad' = pad(p_{N-1}) im System N (Induktionsschritt), 'kneser' = q_N (N ≤ 75).
Gespeichert (float64, npz):  d = J(x0)^{-1} F(x0)  (Newton-Schritt)  und  Z[:, (i,j)] = J(x0)^{-1} F''(x0)(e_i, e_j), i ≤ j.
  α_w = Σ_k |d_k| w_k,   ω_w = max_{i≤j} Σ_k |Z_{k,ij}| w_k / (w_i w_j)  (exakte Bilinearnorm auf ℓ1_w × ℓ1_w),   h_w = α_w ω_w.
F''(e_i, e_j) = [ i p^{i-1} x^j + j p^{j-1} x^i + p''(p) x^{i+j} ]_{≤N}.
Auswertung (feste Gewichtsfamilien und optimierte Gewichte): d03b_gewichte.py.
"""
import json
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [10, 20, 30, 40, 60]
MODUS = sys.argv[2] if len(sys.argv) > 2 else "pad"
BITS = 1200
p = K.loesungen(max(NS) + 1, bits=BITS)
kn = json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_80.json"))
os.makedirs(K.HIER + "/d03_daten", exist_ok=True)


def shift(pol, s, N):
    v = [arb(0)] * (N + 1)
    for k in range(N + 1 - s):
        v[k + s] = pol[k]
    return v


t0 = time.time()
for N in NS:
    ctx.prec = 400 + 3 * N
    if MODUS == "pad":
        a = [arb(x.mid()) for x in p[N - 1]] + [arb(0)]
    elif MODUS == "selbst":
        a = [arb(x.mid()) for x in p[N]]
    else:
        a = [arb(s) for s in kn[:N + 1]]
    P = K.potenzen(a, N)
    F, J = K.FJ(a, N, P)
    g = arb_poly([0])
    for m in range(2, N + 1):
        g += m * (m - 1) * a[m] * P[m - 2]
    paare = [(i, j) for i in range(N + 1) for j in range(i, N + 1)]
    Y = arb_mat(N + 1, len(paare))
    for c, (i, j) in enumerate(paare):
        v = shift(g, i + j, N) if i + j <= N else [arb(0)] * (N + 1)
        if i >= 1:
            u = shift(P[i - 1], j, N)
            v = [v[k] + i * u[k] for k in range(N + 1)]
        if j >= 1:
            u = shift(P[j - 1], i, N)
            v = [v[k] + j * u[k] for k in range(N + 1)]
        for k in range(N + 1):
            Y[k, c] = v[k]
    Z = J.solve(Y, algorithm="approx")
    d = K.loese(J, F)
    if MODUS == "selbst":
        # reduzierter Schritt: Tiefteil-System F_N(a) + t b(a) = 0 mit t = a_{N+1}^{(N+1)}, b = [p_N^{N+1}]_{<=N}
        b = (P[N] * arb_poly(a)).truncate(N + 1)
        t = p[N + 1][N + 1]
        d = K.loese(J, [t * b[k] for k in range(N + 1)])
    Zf = np.array([[float(Z[k, c].mid()) for c in range(len(paare))] for k in range(N + 1)])
    df = np.array([float(x.mid()) for x in d])
    np.savez_compressed(K.HIER + f"/d03_daten/{MODUS}_{N}.npz", Z=Zf, d=df, paare=np.array(paare))
    print(f"N={N} fertig ({time.time()-t0:.0f}s)", flush=True)
