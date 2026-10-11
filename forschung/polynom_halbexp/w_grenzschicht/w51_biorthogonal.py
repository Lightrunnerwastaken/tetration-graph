"""W51 — Ist die Eckmode das LETZTE biorthogonale Paar der Bilinearform ⟨⟨A, Q⟩⟩ = (1/2πi)∮ A(1/x)·[Q(f(x)) + m(x)Q(x)] dx/x (Matrix J_N)?
Q_N (Grad N, vorwärts): exakte Interpolation Q_N∘f + m·Q_N = O(x^N) bei 0, d. h. die ersten N Zeilen von J_N annullieren q (q_N = 1).
A_N (Grad N in 1/x, adjungiert): aᵀ J_N[:, :N] = 0 (a_N = 1).  h_N := aᵀ J_N q.
Vergleich (ρ-skaliert, ℓ²_ρ) mit den Eckmoden u₁ (vorwärts) und v₁ (adjungiert) aus w01; Test σ₁ ≈ ‖q‖‖a‖/|h_N|.
Aufruf: python w51_biorthogonal.py N1,N2,...
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho = W.RHO
for N in [int(v) for v in sys.argv[1].split(",")]:
    t0 = time.time()
    n = N + 1
    ctx.prec = int(1.5 * N) + 300
    a_ = W.kubik_a(N, coef)
    J = W.J_mat(a_, N)
    # Q_N: J[0:N, 0:N] q' = −J[0:N, N]
    A1 = arb_mat(N, N, [J[i, j] for i in range(N) for j in range(N)])
    b1 = arb_mat(N, 1, [-J[i, N] for i in range(N)])
    qp = A1.solve(b1, algorithm="approx")
    q = [qp[i, 0] for i in range(N)] + [arb(1)]
    # A_N: J[0:N, 0:N]ᵀ a' = −J[N, 0:N]ᵀ
    b2 = arb_mat(N, 1, [-J[N, j] for j in range(N)])
    ap = A1.transpose().solve(b2, algorithm="approx")
    a = [ap[i, 0] for i in range(N)] + [arb(1)]
    hN = sum((a[i] * sum((J[i, j] * q[j] for j in range(n)), arb(0)) for i in range(n)), arb(0))
    # ρ-Skalierung: vorwärts δ_k ρ^k, adjungiert y_k ρ^{−k}
    qs = np.array([float((q[k] * arb(rho) ** k).mid()) for k in range(n)])
    as_ = np.array([float((a[k] * arb(rho) ** (-k)).mid()) for k in range(n)])
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, U, Vt = d["S"], d["U"], d["Vt"]
    nq, na = np.linalg.norm(qs), np.linalg.norm(as_)
    ovq = [abs(np.dot(qs / nq, U[:, i])) for i in range(2)]
    ova = [abs(np.dot(as_ / na, Vt[i])) for i in range(2)]
    h = float(hN.mid())
    print(f"N={N} ({time.time()-t0:.0f}s): |⟨Q_N, u₁⟩|, |⟨Q_N, u₂⟩| = {ovq[0]:.4f}, {ovq[1]:.4f}; |⟨A_N, v₁⟩|, |⟨A_N, v₂⟩| = {ova[0]:.4f}, {ova[1]:.4f}; "
          f"‖q‖‖a‖/|h_N| = {nq*na/abs(h):.4g} gegen σ₁ = {S[0]:.4g}", flush=True)
