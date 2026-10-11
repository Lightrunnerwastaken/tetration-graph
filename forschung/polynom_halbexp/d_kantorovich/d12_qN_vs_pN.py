"""D12 — Warum ist J_N(q_N)^{-1} gross, J_N(p_N)^{-1} nicht?  Kleinste Singulärwerte der umskalierten Inversen
(ℓ2 mit w_k = |L|^k (1+k)^s) für q_N (Kneser-Taylor) und p_N, und Vergleich des Singulärvektors mit den Störrichtungen δ_{±1}.
"""
import json
import sys

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [100, 150, 194, 250]
kn = json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_1000.json"))["a"]
dd = json.load(open(K.GEMEINSAM + "/stoerrichtungen_delta.json"))["delta"]["1"]
d1 = np.array([float(x) for x in dd["re"]]) + 1j * np.array([float(x) for x in dd["im"]])
p = K.loesungen(max(NS), bits=1200)
L = 1.3745570107437075


def inv_resc(a, N):
    ctx.prec = 500 + 3 * N
    _, J = K.FJ(a, N)
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    return np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]


for N in NS:
    k = np.arange(N + 1)
    for name, a in (("q_N", [arb(str(x)) for x in kn[:N + 1]]), ("p_N", [arb(x.mid()) for x in p[N]])):
        G = inv_resc(a, N)
        for s in (0, 4):
            w = (1.0 + k) ** s
            A = G * w[:, None] / w[None, :]
            U, S, Vt = np.linalg.svd(A)
            u = U[:, 0] / w           # grösster Singulärvektor im Bild (umskaliert, ungewichtet)
            # Vergleich mit δ_{±1}: Anpassung von u an Re/Im von T_N δ_1 (umskaliert)
            M = np.vstack([np.real(d1[:N + 1]) * L ** k, np.imag(d1[:N + 1]) * L ** k]).T
            c, *_ = np.linalg.lstsq(M, u, rcond=None)
            rest = np.linalg.norm(u - M @ c) / np.linalg.norm(u)
            print(f"N={N} {name} s={s}: σ_max(J⁻¹) = {S[0]:.4g}, nächste {S[1]:.4g}, {S[2]:.4g};  Singulärvektor vs δ±1: Rest {rest:.3f}",
                  flush=True)
