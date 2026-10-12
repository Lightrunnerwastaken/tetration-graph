"""k16 — Fast-Orthogonalität der Köpfe (Cotlar–Stein). T_n = c_n ⊗ r̂_n (w-Kopfterme, Ĵ_N⁻¹ = Σ_{n≤N} T_n).
Cotlar–Stein: ‖Σ T_n‖ ≤ max(sup_n Σ_m ‖T_n*T_m‖^{1/2}, sup_n Σ_m ‖T_nT_m*‖^{1/2}),
‖T_n*T_m‖ = |⟨c_n,c_m⟩|·‖r̂_n‖‖r̂_m‖, ‖T_nT_m*‖ = ‖c_n‖‖c_m‖·|⟨r̂_n,r̂_m⟩|.
Gemessen: beide Cotlar–Stein-Summen für n ≤ N (N = 100, 200, 300) und die Abklingform der Gram-Kosinusse in |n − m|.
Daten: k13_koepfe_2d_kubisch_300.npz (Köpfe aus k13)."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "k13_koepfe_2d_kubisch_300.npz"); C = d["c"]; R = d["r"]; Nr = C.shape[0] - 1
nc = np.linalg.norm(C, axis=1); nr = np.linalg.norm(R, axis=1)
Gc = C @ C.T; Gr = R @ R.T
cosc = Gc / np.outer(nc, nc); cosr = Gr / np.outer(nr, nr)
for N in sorted(set([100, 200, 300, 400, 500, Nr])):
    A1 = np.sqrt(np.abs(Gc[:N + 1, :N + 1]) * np.outer(nr[:N + 1], nr[:N + 1])).sum(axis=1)
    A2 = np.sqrt(np.outer(nc[:N + 1], nc[:N + 1]) * np.abs(Gr[:N + 1, :N + 1])).sum(axis=1)
    S = sum(np.outer(C[n], R[n]) for n in range(N + 1))
    print("N = %d: Cotlar–Stein-Summen sup_n: %.4g (c-Gram), %.4g (r-Gram); tatsächlich ‖Ĵ_N⁻¹‖ = %.6f; Σ‖T_n‖ = %.4g"
          % (N, A1.max(), A2.max(), np.linalg.norm(S, 2), (nc[:N + 1] * nr[:N + 1]).sum()))
for n in [150, 250, Nr - 50]:
    ks = [1, 2, 3, 5, 10, 20, 40, 80]
    print("n = %d: |cos(c_n, c_{n−k})|: %s" % (n, "  ".join("k=%d: %.3f" % (k, abs(cosc[n, n - k])) for k in ks)))
    print("        |cos(r̂_n, r̂_{n−k})|: %s" % ("  ".join("k=%d: %.3f" % (k, abs(cosr[n, n - k])) for k in ks)))
