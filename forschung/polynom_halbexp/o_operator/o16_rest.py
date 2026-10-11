"""O16 — Woher kommt die Restkonstante 1 + δ (radial 1.046e-3, kubisch 5.3e-5)? ‖P_n Y_N P_n‖ (Niederfrequenzblock) gegen σ_{r+1}(Y_N),
und Lage des Rest-Singulärvektors (k-Verteilung, Ort auf dem Kreis)."""
import sys
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]; N = int(sys.argv[2])
Y = K.Y_matrix(name, N, bits=int(0.85 * N) + 200)
U, s, Vh = np.linalg.svd(Y)
r = int((s > 1.01).sum())
print(f"{name} N={N}: σ_1..r+3 = {np.array2string(s[:r+3], precision=7)}")
for n in (N // 8, N // 4, N // 2, 3 * N // 4):
    sb = np.linalg.svd(Y[:n + 1, :n + 1], compute_uv=False)
    print(f"  ‖P_n Y P_n‖ (n = {n}): {sb[0]:.7f};  ‖Y P_n‖ = {np.linalg.svd(Y[:, :n+1], compute_uv=False)[0]:.7f}")
k = np.arange(N + 1)
phi = np.linspace(-np.pi, np.pi, 4001)
for i in (r, r + 1):
    for nm, vec in (("Eingang v", Vh[i].conj()), ("Ausgang u", U[:, i])):
        w = np.abs(vec) ** 2
        vals = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * phi), vec))
        print(f"  Rest-Vektor i={i+1} (σ−1 = {s[i]-1:.3e}) {nm}: k-Schwerpunkt {(w*k).sum():.1f}, Masse k<20: {w[:20].sum():.3f}, "
              f"Max auf Kreis bei {np.degrees(phi[np.argmax(vals)]):.1f}°")
