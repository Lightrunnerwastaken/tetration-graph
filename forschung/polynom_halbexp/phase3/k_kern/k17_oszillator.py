"""k17 — Woher kommt die Fast-Orthogonalität der Köpfe? Test der Oszillator-Struktur T_{n+1} − 2cos(ω)T_n + T_{n−1} ≈ 0.
T_n = c_n ⊗ r̂_n (w-Kopfterme, Daten k13_koepfe_2d_kubisch_300.npz). Frobenius-Norm über Gram-Matrizen (exakt für Summen von Rang-1-Termen).
Gemessen: optimales ω je Fenster, relativer Defekt ε_n = ‖T_{n+1} − 2cosω T_n + T_{n−1}‖/‖T_n‖, Skalierung von ε_n mit n; Vergleich mit zwei Frequenzen."""
import sys
import numpy as np
from scipy.optimize import minimize_scalar
sys.stdout.reconfigure(encoding="utf-8")
d = np.load("k13_koepfe_2d_kubisch_300.npz"); C = d["c"]; R = d["r"]; Nr = C.shape[0] - 1
Gc = C @ C.T; Gr = R @ R.T
def fro2(coef, idx):
    s = 0.0
    for a, i in zip(coef, idx):
        for b, j in zip(coef, idx):
            s += a * b * Gc[i, j] * Gr[i, j]
    return s
def eps(n, w):
    return np.sqrt(max(fro2([1, -2 * np.cos(w), 1], [n + 1, n, n - 1]), 0)) / np.sqrt(fro2([1], [n]))
for lo, hi in [(100, 150), (150, 200), (200, 250), (250, 299)]:
    f = lambda w: sum(eps(n, w) ** 2 for n in range(lo, hi))
    r = minimize_scalar(f, bounds=(0.01, np.pi - 0.01), method="bounded")
    w = r.x
    e = [eps(n, w) for n in range(lo, hi)]
    print("Fenster %d–%d: ω = ±%.5f (mod 2π), mittlerer Defekt ε = %.4f, n·ε = %.3f; ohne Oszillator (ω=0, zweite Differenz): %.4f"
          % (lo, hi, w, np.mean(e), np.mean([n * x for n, x in zip(range(lo, hi), e)]), np.mean([eps(n, 0.0) for n in range(lo, hi)])))
# zwei Frequenzen: (S − e^{iω1})(S − e^{−iω1})(S − e^{iω2})(S − e^{−iω2}) T_n — reelle Koeffizienten
def eps4(n, w1, w2):
    p = np.polymul([1, -2 * np.cos(w1), 1], [1, -2 * np.cos(w2), 1])   # Koeffizienten für T_{n+2} … T_{n−2}
    return np.sqrt(max(fro2(list(p), [n + 2, n + 1, n, n - 1, n - 2]), 0)) / np.sqrt(fro2([1], [n]))
for lo, hi in [(150, 200), (250, 297)]:
    best = None
    for w1 in np.linspace(1.3, 1.8, 26):
        for w2 in np.linspace(0.05, 3.1, 62):
            v = sum(eps4(n, w1, w2) ** 2 for n in range(lo, hi, 3))
            if best is None or v < best[0]:
                best = (v, w1, w2)
    v, w1, w2 = best
    e = [eps4(n, w1, w2) for n in range(lo, hi)]
    print("zwei Frequenzen, Fenster %d–%d: ω1 = %.3f, ω2 = %.3f, mittlerer Defekt %.4f" % (lo, hi, w1, w2, np.mean(e)))
