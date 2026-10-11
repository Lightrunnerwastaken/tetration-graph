"""O37 — Vergleich der Innenschicht-Profile aus o36 über N: Koeffizientenprofil |v_k| in κ = k/N und in ℓ = N − k,
Kreisprofil |v| in y = N(θ_L − φ), Überlaufprofil log10|[J∞v]_{N+j}| in j und j/N. Aufruf: python o37_profilvergleich.py modell N1,N2,...
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
Ns = [int(x) for x in sys.argv[2].split(",")]
D = {N: np.load(f"o36_{name}_{N}.npz") for N in Ns}
for mode in range(2):
    print(f"===== Mode {mode}")
    print("s:", " ".join(f"N={N}:{D[N]['s'][mode]:.3f}" for N in Ns))
    print("Kreisprofil |v(y)|/max, y = N(θ_L−φ):")
    for y in (-200, -100, -50, -20, -10, -5, 0, 5, 10, 15, 20, 30, 40):
        row = []
        for N in Ns:
            ys, pr = D[N]["ys"], D[N]["prof"][mode]
            row.append(pr[np.argmin(np.abs(ys - y))] / pr.max())
        print(f"  y = {y:+5d}: " + " ".join(f"{v:.3e}" for v in row))
    print("Koeffizienten |v_k| (unit) bei ℓ = N − k:")
    for l in (0, 1, 2, 3, 5, 10, 20, 40, 80):
        print(f"  ℓ = {l:3d}: " + " ".join(f"{abs(D[N]['V'][mode][N - l]):.3e}" for N in Ns))
    print("Koeffizienten |v_k| bei κ = k/N:")
    for kap in (0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 1.0):
        print(f"  κ = {kap:.2f}: " + " ".join(f"{abs(D[N]['V'][mode][int(round(kap * N))]):.3e}" for N in Ns))
    print("Überlauf log10|[J∞v]_(N+j)| bei j:")
    for j in (-3, -2, -1, 0, 1, 2, 3, 5, 10, 20):
        print(f"  j = {j:+3d}: " + " ".join(f"{np.log10(abs(D[N]['JV'][mode][N + j]) + 1e-300):8.2f}" for N in Ns))
    print("Überlauf bei j/N:")
    for jn in (0.01, 0.02, 0.05, 0.1, 0.2):
        print(f"  j/N = {jn:.2f}: " + " ".join(f"{np.log10(abs(D[N]['JV'][mode][N + int(round(jn * N))]) + 1e-300) / N:8.4f}" for N in Ns) + "   (log10/N)")
