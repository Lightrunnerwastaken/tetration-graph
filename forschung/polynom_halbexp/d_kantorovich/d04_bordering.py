"""D04 — Gleichmässige Stabilität |J_N^{-1}|_w über die Bordering-Rekursion.

Gewichte w_k = |L|^k (1+k)^s (gleich in Urbild und Bild), Operatornorm auf ℓ1_w = max_j |J^{-1} e_j|_w / w_j.
Zerlegung des Schritts N -> N+1 in zwei Teile:
  (A) Rand: J_{N+1}(pad p_N) = [[J_N, b],[cᵀ, d]]  =>  J_{N+1}(pad)^{-1} = [[J_N^{-1} + u vᵀ/σ, −u/σ], [−vᵀ/σ, 1/σ]],
      u = J_N^{-1} b, vᵀ = cᵀ J_N^{-1}, σ = d − cᵀu (Schur-Komplement).
      Schranke: |J_{N+1}(pad)^{-1}| ≤ |J_N^{-1}| + ε_N,  ε_N = |v|_{w*} (|u|_w + w_{N+1}) / |σ|,  |v|_{w*} = max_j |v_j|/w_j
      (und für die letzte Spalte (|u|_w + w_{N+1})/(|σ| w_{N+1})).
  (B) Störung: J_{N+1}(p_{N+1}) = J_{N+1}(pad) + ΔJ;  π_N = |J_{N+1}(pad)^{-1} ΔJ|_w  =>  |J(p_{N+1})^{-1}| ≤ |J(pad)^{-1}|/(1 − π_N).
Ausgabe: d04_bordering.json.
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
SS = [0, 2, 4, 5]
p = K.loesungen(NMAX + 1, bits=1200)
absL = 1.3745570107437075


def inv_float(J):
    n = J.nrows()
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    return np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]), Ji


def gew(n, s):
    k = np.arange(n)
    return absL ** k * (1.0 + k) ** s


def opnorm(M, w):
    """ℓ1_w-Operatornorm: max_j Σ_i |M_ij| w_i / w_j  (und argmax j)."""
    col = (np.abs(M) * w[:, None]).sum(axis=0) / w
    return float(col.max()), int(col.argmax())


aus = {}
t0 = time.time()
Jinv_cache = None
for N in range(4, NMAX + 1):
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    a1 = [arb(x.mid()) for x in p[N + 1]]
    ap = a + [arb(0)]
    _, JN = K.FJ(a, N)
    _, Jpad = K.FJ(ap, N + 1)
    _, Jnew = K.FJ(a1, N + 1)
    JNi, _ = inv_float(JN)
    Jpi, Jpi_arb = inv_float(Jpad)
    Jni, _ = inv_float(Jnew)
    # Randgrössen
    b = np.array([float(Jpad[k, N + 1].mid()) for k in range(N + 1)])
    c = np.array([float(Jpad[N + 1, j].mid()) for j in range(N + 1)])
    # u, v aus der exakten Blockinversen (stabiler als float-Produkte): letzte Spalte/Zeile von Jpad^{-1}
    sigma = 1.0 / Jpi[N + 1, N + 1]
    u = -Jpi[:N + 1, N + 1] * sigma
    v = -Jpi[N + 1, :N + 1] * sigma
    # Störung
    dJ = Jnew - Jpad
    P = Jpi_arb * dJ
    Pf = np.array([[float(P[i, j].mid()) for j in range(N + 2)] for i in range(N + 2)])
    zeile = {"sigma": sigma}
    for s in SS:
        wN, wN1 = gew(N + 1, s), gew(N + 2, s)
        nN, jN = opnorm(JNi, wN)
        npad, jpad = opnorm(Jpi, wN1)
        nnew, jnew = opnorm(Jni, wN1)
        vstar = float(np.max(np.abs(v) / wN))
        uw = float((np.abs(u) * wN).sum())
        eps = vstar * (uw + wN1[-1]) / abs(sigma)
        letzte = (uw + wN1[-1]) / (abs(sigma) * wN1[-1])
        pi, _ = opnorm(Pf, wN1)
        zeile[str(s)] = {"JN": nN, "jN": jN, "Jpad": npad, "Jnew": nnew, "jnew": jnew, "eps": eps, "vstar": vstar,
                         "uw": uw, "letzte": letzte, "pi": pi}
    aus[N] = zeile
    if N % 10 == 0 or N < 8:
        print(f"N={N:3d} σ={sigma:.3e} ({time.time()-t0:.0f}s)")
        for s in SS:
            z = zeile[str(s)]
            print(f"   s={s}: |J_N^-1|={z['JN']:.4f} (j={z['jN']})  |J_pad^-1|={z['Jpad']:.4f}  |J_N+1^-1|={z['Jnew']:.4f}"
                  f"  ε={z['eps']:.3e} [v*={z['vstar']:.2e} |u|={z['uw']:.2e}]  letzte={z['letzte']:.2e}  π={z['pi']:.2e}",
                  flush=True)
json.dump({str(k): v for k, v in aus.items()}, open(K.HIER + "/d04_bordering.json", "w"))
