"""O17 — Restkonstante δ₀ = ‖G‖ − 1 für Kneser (A_f aus den Koeffizienten in gemeinsam/, kein Motor): ‖P_n Y_N P_n‖ für kleine n, N = 60…140,
dazu σ₁..σ₄(Y_N). Präzision: Koeffizientenfehler 4e-47 ⇒ N ≲ 150."""
import sys
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
for N in [int(s) for s in sys.argv[1].split(",")]:
    Y = K.Y_matrix("kneser", N, bits=400)
    s = np.linalg.svd(Y, compute_uv=False)
    bl = [np.linalg.svd(Y[:n + 1, :n + 1], compute_uv=False)[0] - 1 for n in (20, 30, 40, 50)]
    print(f"kneser N={N}: σ−1 = {np.array2string(s[:5]-1, precision=4)}; ‖P_nYP_n‖−1 (n=20/30/40/50) = {np.array2string(np.array(bl), precision=6)}", flush=True)
