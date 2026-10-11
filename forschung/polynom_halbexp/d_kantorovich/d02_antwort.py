"""D02 — Antwort der tiefen Koeffizienten auf die hohen: R_N = J_N^{-1} B_N.

Exakt (Lemma D1): Für p = p_lo + p_hi (p_hi = O(x^{N+1})) ist
    [p∘p − exp]_{≤N} = F_N(a_lo) + B_N(a_lo) a_hi,   B_N(a)_{j,k} = [x^j] p_lo^k  (k > N).
Also ist Knesers q_N die exakte Lösung des N-Systems mit „Schwanzkraft“ B_N q_hi, und p_N − q_N ≈ J_N^{-1} B_N q_hi.
Gemessen: ρ_k^{(N)} = |J_N^{-1} b_k|_X / |L|^k  (X = Σ|v_j||L|^j bzw. sup auf |x| = |L|), k = N+1 … KMAX·N.
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
BITS = 1200
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [10, 20, 30, 40, 60, 80, 100, 120, 150]
FAKTOR = 4
p = K.loesungen(max(NS), bits=BITS)
ctx.prec = BITS
absL = abs(K.L_wert())
absLf = float(absL.mid())


def sup_skaliert(v_scaled):
    """v_scaled[j] = v_j |L|^j (float) -> max_phi |Σ v_scaled[j] e^{ijφ}|."""
    n = 2048
    return float(np.max(np.abs(np.fft.fft(np.asarray(v_scaled, dtype=float), n))))


aus = {}
t0 = time.time()
for N in NS:
    a = p[N]
    KM = FAKTOR * N + 10
    P = K.potenzen(a, KM)          # p^k mod x^{KM+1}; wir brauchen nur [x^j], j ≤ N
    F, J = K.FJ(a, N)
    ks = list(range(N + 1, KM + 1))
    B = arb_mat(N + 1, len(ks))
    for c, k in enumerate(ks):
        for j in range(N + 1):
            B[j, c] = P[k][j]
    X = J.solve(B, algorithm="approx")
    Lpot = [absL ** j for j in range(KM + 2)]
    zeilen = []
    for c, k in enumerate(ks):
        vs = [X[j, 0 + c] * Lpot[j] for j in range(N + 1)]
        l1 = sum((abs(v) for v in vs), arb(0)) / Lpot[k]
        sup = sup_skaliert([float((v / Lpot[k]).mid()) for v in vs])
        bnorm = sum((abs(B[j, c]) * arb("0.5") ** j for j in range(N + 1)), arb(0))
        zeilen.append((k, float(l1.mid()), sup, float(bnorm.mid())))
    aus[N] = zeilen
    kk = [z[0] for z in zeilen]
    l1s = np.array([z[1] for z in zeilen])
    sups = np.array([z[2] for z in zeilen])
    i1 = 0
    print(f"N={N:3d}  rho_(N+1): l1={l1s[0]:.3e} sup={sups[0]:.3e} | max_k l1={l1s.max():.3e} bei k={kk[int(l1s.argmax())]} "
          f"| rho_(2N)={l1s[N-1]:.3e} | rho_(4N)={l1s[min(len(l1s)-1, 3*N-1)]:.3e}  ({time.time()-t0:.0f}s)", flush=True)
json.dump({str(N): v for N, v in aus.items()}, open(K.HIER + "/d02_antwort.json", "w"))
