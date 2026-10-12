"""Strang V — (S_γ) für das ECHTE p_N statt f_K (Auftrag Leitung 22:52, Punkt 2): gilt sup_N ‖(Ĵ_N^{p_N})⁻¹‖ < ∞ ?

Ĵ_N^{P} := T_N^w[(δ∘P + m_P·δ)∘k] mit der konformen Abbildung k aus v12 (Linse S, Ks Profilform), Spalten wie v12.
  P = p_N (Strang A, gemeinsam/pN_N100_bis_1000.json, 60 Stellen, Kneser-Zweig), m_P = p_N'∘p_N (Tangente des Newton-Verfahrens,
      exakt als Polynomkomposition modulo x^{N+1});
  Vergleich P = f_K (gemeinsam/kneser_halbexp_koeff_4000.json, 50 Stellen), m = f_K'∘f_K aus m·f_K' = e^x (wie v12 kneser_z).
Hinweis: T_N^w[g∘k] hängt nur von T_N^x g ab (k(0) = 0) — der Vergleich misst also genau η = p_N − T_N f_K und m_P − T_N m.
Ausgabe je N: ‖(Ĵ^{p_N})⁻¹‖, ‖(Ĵ^{f_K})⁻¹‖, Δ_inv = ‖(Ĵ^{p_N})⁻¹ − (Ĵ^{f_K})⁻¹‖ und ε_N = Σ_{k≤N}|p_k − a_k|ρ^k (davon Band k ≥ 0.8N).
Aufruf: python v20_pN_S.py N1,N2,... S
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

ARGS = sys.argv[1:]
Ns = [int(x) for x in ARGS[0].split(",")]
S = ARGS[1]
sys.argv = ["v12", "kneser_z", "100", S]
import v12_konform_S as V12  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
GEM = Path(__file__).resolve().parents[2] / "gemeinsam"
PN = json.loads((GEM / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"]
AK = json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]


def koeff_pN(N):
    p = [arb(x) for x in PN[str(N)]]
    assert len(p) == N + 1
    P = arb_poly(p)
    dP = P.derivative()
    r = arb_poly([0])
    for c in reversed(dP.coeffs()):                       # p'∘p mod x^{N+1} (Horner, exakte Abschneidung)
        r = (r * P).truncate(N + 1) + c
    m = [r[k] if k <= r.degree() else arb(0) for k in range(N + 1)]
    return p, m


def koeff_fK(N):
    a = [arb(AK[k]) for k in range(N + 2)]
    fp = [a[k + 1] * (k + 1) for k in range(N + 1)]
    ex = [arb(1)]
    for k in range(1, N + 1):
        ex.append(ex[-1] / k)
    m = [arb(0)] * (N + 1)
    for k in range(N + 1):
        s_ = ex[k]
        for j in range(k):
            s_ -= m[j] * fp[k - j]
        m[k] = s_ / fp[0]
    return a[:N + 1], m


def inverse(J):
    n = J.nrows()
    Id = arb_mat(n, n)
    for i in range(n):
        Id[i, i] = 1
    return J.solve(Id, algorithm="approx")


def zu_np(M):
    return np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


def main():
    F, it, d = V12.theodorsen()
    print(f"Kneser, Linse S = {S} (Ks Profilform), Theodorsen {it} It.", flush=True)
    for N in Ns:
        t0 = time.time()
        bits = int(1.2 * N) + 400
        ctx.prec = bits
        rho = abs(V12.L_arb())
        kp, kf = koeff_pN(N), koeff_fK(N)
        eps = [float((abs(kp[0][k] - kf[0][k]) * rho ** k).mid()) for k in range(N + 1)]
        out = {}
        for name, kk in (("pN", kp), ("fK", kf)):
            V12.modell_koeff = lambda n_, kk=kk: kk
            J = V12.J_hat(N, F, bits)
            ctx.prec = bits
            out[name] = inverse(J)
        Xp, Xf = zu_np(out["pN"]), zu_np(out["fK"])
        ctx.prec = bits
        dX = zu_np(out["pN"] - out["fK"])
        sp, sf = np.linalg.norm(Xp, 2), np.linalg.norm(Xf, 2)
        dn = np.linalg.norm(dX, 2)
        e_all, e_top = sum(eps), sum(eps[int(0.8 * N):])
        print(f"N = {N:5d}  ‖(Ĵ^pN)⁻¹‖ = {sp:.6f}  ‖(Ĵ^fK)⁻¹‖ = {sf:.6f}  Δ_inv = {dn:.3e}  ε_N = {e_all:.3e} (Band ≥ 0.8N: {e_top:.3e})  "
              f"Δ_inv/ε_N = {dn / e_all:.3e}  [{bits} Bit, {time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
