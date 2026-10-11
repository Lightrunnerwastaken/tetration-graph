"""O05 — Ist K_N (ℓ²_ρ) akkretiv, Re⟨K_N δ, δ⟩ ≥ 0? (würde ‖(I+K_N)^{-1}‖ ≤ 1 sofort geben). Modelle elliptisch/kubisch."""
import sys
import numpy as np
from flint import arb_mat, acb_mat
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]; mu = float(sys.argv[3]) if len(sys.argv) > 3 else None
for N in [int(s) for s in sys.argv[2].split(",")]:
    J, M, kp = K.matrizen(name, N, bits=int(1.3*N) + 250, mu=mu)
    Km = M.solve(J, algorithm="approx")          # M^{-1} J = I + K
    A = K.to_np(Km, kp)
    Kn = A - np.eye(N + 1)
    H = (Kn + Kn.conj().T) / 2
    ev = np.linalg.eigvalsh(H)
    G = A.conj().T @ A - np.eye(N + 1)
    eg = np.linalg.eigvalsh(G)
    print(f"{name} {mu} N={N}: Re K: min {ev[0]:+.3e}, #neg {int((ev < -1e-12*abs(ev).max()).sum())}, max {ev[-1]:.2e} | A*A−I: min {eg[0]:+.3e}", flush=True)
