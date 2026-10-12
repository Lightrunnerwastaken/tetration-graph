"""k20 — Oszillator-Struktur NUR im Eckanteil: T^E_n = c_n ⊗ r̂^E_n, r̂^E_n = adjungierter w-Kopf, auf ein festes Winkelfenster um ±θ_c beschränkt.
Defekt ε_n(ω) = ‖T^E_{n+1} − 2cos ω T^E_n + T^E_{n−1}‖/‖T^E_n‖ mit ω = Φ(θ_c) (Phaseninkrement an der Ecke, k19) und mit optimalem ω.
Klein und ~1/n ⇒ Summation durch Teile im Eckanteil trägt (log N). Daten: k13_koepfe_2d_kubisch_300.npz."""
import sys
import numpy as np
from scipy.optimize import minimize_scalar
sys.stdout.reconfigure(encoding="utf-8")
d = np.load("k13_koepfe_2d_kubisch_300.npz"); C = d["c"]; R = d["r"]; Nr = R.shape[0] - 1
M = 4096
th = 2 * np.pi * np.arange(M) / M - np.pi
E = np.exp(-1j * np.outer(np.arange(Nr + 1), th))
G = R @ E
z = sum(G[n] * np.conj(G[n - 1]) for n in range(270, 301)); Phi = np.angle(z)
pk = th[np.argmax(np.abs(G[300]) * (th > 0))]
ic = np.argmin(np.abs(th - pk)); omega_c = Phi[ic]
print("Ecke bei θ_c = %.2f°, Phaseninkrement ω_c = Φ(θ_c) = %.5f rad (Gs Wert an L: −1.48542)" % (np.degrees(pk), omega_c))
for half in [3.0, 8.0]:
    chi = ((np.abs(th - pk) <= np.radians(half)) | (np.abs(th + pk) <= np.radians(half))).astype(float)
    # zurück in Koeffizienten (Länge Nr+1, Projektion auf die Indizes ≤ n)
    RE = np.zeros_like(R, dtype=complex)
    for n in range(Nr + 1):
        coef = (G[n] * chi) @ np.exp(1j * np.outer(th, np.arange(Nr + 1))) / M
        coef[n + 1:] = 0
        RE[n] = coef
    Gc = C @ C.T; Gr = RE @ RE.conj().T
    def fro(coef, idx):
        s = 0
        for a, i in zip(coef, idx):
            for b, j in zip(coef, idx):
                s += a * np.conj(b) * Gc[i, j] * Gr[j, i]
        return np.sqrt(abs(s))
    def eps(n, w): return fro([1, -2 * np.cos(w), 1], [n + 1, n, n - 1]) / fro([1], [n])
    frac = np.mean([np.linalg.norm(RE[n]) ** 2 / np.linalg.norm(R[n]) ** 2 for n in range(250, 300)])
    for lo, hi in [(150, 200), (250, 299)]:
        e_c = [eps(n, omega_c) for n in range(lo, hi)]
        r = minimize_scalar(lambda w: sum(eps(n, w) ** 2 for n in range(lo, hi)), bounds=(0.01, np.pi - 0.01), method="bounded")
        e_o = [eps(n, r.x) for n in range(lo, hi)]
        print("Fenster ±%.0f° (Massenanteil %.3f), n = %d–%d: Defekt mit ω_c: %.4f (n·ε = %.2f); optimales ω = ±%.4f: %.4f (n·ε = %.2f)"
              % (half, frac, lo, hi, np.mean(e_c), np.mean([n * x for n, x in zip(range(lo, hi), e_c)]), r.x, np.mean(e_o), np.mean([n * x for n, x in zip(range(lo, hi), e_o)])))
    # Summation durch Teile direkt: ‖Σ_{n≤N} T^E_n‖ vs Σ‖T^E_n‖
    for N in [150, 200, 250, 300]:
        S = sum(np.outer(C[n], RE[n]) for n in range(N + 1))
        print("  N = %d: ‖Σ T^E_n‖ = %.4f, Σ‖T^E_n‖ = %.3f" % (N, np.linalg.norm(S, 2), sum(np.linalg.norm(C[n]) * np.linalg.norm(RE[n]) for n in range(N + 1))))
