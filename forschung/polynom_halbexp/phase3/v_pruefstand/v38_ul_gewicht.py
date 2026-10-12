"""Strang V — VE4 bis N = 2000 (Auftrag Leitung 02:05): UL-Zerlegung der Inversen X = Ĵ_N⁻¹ = U'L' (U' oben, L' unten mit Einsen) per
Rückwärts-Elimination auf X selbst (Satz K5: X − Kopf_N = Ĵ_{N−1}⁻¹ ⊕ 0, also ist Spalte k von U' = Ĵ_k⁻¹e_k, Zeile k von L' = e_kᵀĴ_k⁻¹/(Ĵ_k⁻¹)_kk).
Arithmetik: X aus arb_mat.solve (1.2N + 400 Bit), Elimination in Mittelpunkt-Arithmetik mit P Bit (Dynamik der Köpfe ~e^{±0.07N}), Rang-1-Schritte
als C-Matrixoperationen, alle 250 Schritte auf den führenden Block verkleinert.
Dann Diagonalgewicht: s_k = (‖u_k‖/‖ℓ_k‖)^{1/2} (Spalten und Zeilen gleich gross), danach Feinausgleich in Gleitkomma.
Ausgabe: ‖X‖, min_W ‖U'W⁻¹‖·‖WL'‖, Schur-Summen, Profil log s_k (Steigung), ‖u_k‖·‖ℓ_k‖ (Kopfgrössen), Datei mit s_k.
Aufruf: python v38_ul_gewicht.py MODELL N S [tau1] [p]       MODELL kubisch | kneser
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, arb_mat, ctx

ARGS = sys.argv[1:]
MODELL, N, S = ARGS[0], int(ARGS[1]), ARGS[2]
TAU1 = ARGS[3] if len(ARGS) > 3 else "12"
PEXP = ARGS[4] if len(ARGS) > 4 else "1.3"
ITER = ARGS[5] if len(ARGS) > 5 else "300"
sys.argv = ["v20", "100", S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12
V12.S, V12.TAU1, V12.PEXP = float(S), math.radians(float(TAU1)), float(PEXP)
V12.MODELL = "kubisch" if MODELL == "kubisch" else "kneser_z"
sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent


def main():
    t0 = time.time()
    F, it, d = V12.theodorsen()
    bits = int(1.2 * N) + 400
    ctx.prec = bits
    if MODELL == "kneser":
        kk = V20.koeff_fK(N)
        V12.modell_koeff = lambda n_: kk
    J = V12.J_hat(N, F, bits)
    ctx.prec = bits
    X = V20.inverse(J)
    n = N + 1
    t1 = time.time()
    P = int(__import__("os").environ.get("V38_P", int(0.25 * N) + 300))
    ctx.prec = P
    X = X.mid()
    U = np.zeros((n, n))
    L = np.zeros((n, n))
    su = np.zeros(n)
    sl = np.zeros(n)
    m = n                                     # aktuelle Blockgrösse
    for k in range(n - 1, -1, -1):
        if m - (k + 1) >= 250:
            Y = arb_mat(k + 1, k + 1)
            for i in range(k + 1):
                for j in range(k + 1):
                    Y[i, j] = X[i, j]
            X = Y
            m = k + 1
        piv = X[k, k]
        u = arb_mat(m, 1)
        lr = arb_mat(1, m)
        for i in range(k + 1):
            u[i, 0] = X[i, k]
            lr[0, i] = X[k, i] / piv
        # Normen in arb (Dynamik!), dann Gewicht s_k = (‖u‖/‖ℓ‖)^{1/2}
        nu2 = arb(0)
        nl2 = arb(0)
        for i in range(k + 1):
            nu2 += u[i, 0] ** 2
            nl2 += lr[0, i] ** 2
        nu_, nl_ = nu2.sqrt(), nl2.sqrt()
        s = (nu_ / nl_).sqrt()
        su[k] = float(nu_.log().mid())
        sl[k] = float(nl_.log().mid())
        for i in range(k + 1):
            U[i, k] = float((u[i, 0] / s).mid())
            L[k, i] = float((lr[0, i] * s).mid())
        X = (X - u * lr).mid()
    t2 = time.time()
    Xf = U @ L
    nX = np.linalg.norm(Xf, 2)
    a0, b0 = np.linalg.norm(U, 2), np.linalg.norm(L, 2)
    np.savez_compressed(HIER / f"v38_UL_{MODELL}_N{N}_S{S}_t{TAU1}_p{PEXP}.npz", U=U, L=L)
    best = (a0 * b0, a0, b0, np.ones(n))      # Feinoptimierung getrennt: v42_gewicht_opt.py
    prod, a, b, w = best
    A_ = np.abs(U / w[None, :])
    B_ = np.abs(L * w[:, None])
    schurA = math.sqrt(A_.sum(axis=0).max() * A_.sum(axis=1).max())
    schurB = math.sqrt(B_.sum(axis=0).max() * B_.sum(axis=1).max())
    logs = 0.5 * (su - sl) - np.log(w)           # log des Gesamtgewichts W_k (bis auf Konstante)
    kk_ = np.arange(n)
    sel = (kk_ > 0.2 * N) & (kk_ < 0.8 * N)
    steig = np.polyfit(kk_[sel], logs[sel], 1)[0]
    kopf = np.exp(su + sl)                       # ‖u_k‖·‖ℓ_k‖ = Grösse des Kopfterms
    print(f"[{MODELL}, N = {N}, Linse S = {S}, τ₁ = {TAU1}°, p = {PEXP}] ‖Ĵ⁻¹‖ = {nX:.5f}; ‖U'W⁻¹‖·‖WL'‖: ausgeglichen {a0 * b0:.4f}, "
          f"fein {prod:.4f} ({a:.3f} · {b:.3f}); Schur-Summen {schurA:.3f} / {schurB:.3f}; log W_k-Steigung (0.2N…0.8N) = {steig:.5f}; "
          f"Kopfgrösse max {kopf.max():.3e}, Mittel {kopf.mean():.3e}  [Inverse {t1 - t0:.0f} s, UL {t2 - t1:.0f} s, {P} Bit]", flush=True)
    (HIER / f"v38_{MODELL}_N{N}_S{S}_t{TAU1}_p{PEXP}.json").write_text(json.dumps(
        dict(N=N, S=S, logW=[float(x) for x in logs], log_u=[float(x) for x in su], log_l=[float(x) for x in sl], w_fein=[float(x) for x in w],
             prod=prod, a=a, b=b, schur=[schurA, schurB])), encoding="utf-8")


if __name__ == "__main__":
    main()
