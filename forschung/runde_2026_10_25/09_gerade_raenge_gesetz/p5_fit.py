"""Exponent alpha_r aus den theta-Koeffizienten mit Korrekturgliedern (kein Motor).
Modell (aus der lokalen Struktur theta = -z + Theta(u), u = (z-x*)^alpha g(z)):
  a_n e^{2 pi i n x*} n^{1+alpha} = c1 + c2 n^{-alpha} + c3 n^{-1} + c4 n^{-2 alpha}
Rang 8: nur |a_n| aus runde_2026_10_18/02/daten/theta8_20_0.03_64.json (gelesen) -> reelles log-Modell.
Rang 10: komplexe a_n aus daten/theta10_20_0.03_64.json -> komplexes Modell."""
import json
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent
lam = {2: 2.0, 4: 5.171797263545521, 6: 9.49308711431589, 8: 14.7587699373042601}
vor = {8: np.log(lam[4]) / np.log(lam[6]), 10: np.log(lam[6]) / np.log(lam[8]), 6: np.log(2) / np.log(lam[4])}
print("Vorhersagen alpha_r = ln lam_{r-4}/ln lam_{r-2}:", {k: round(v, 9) for k, v in vor.items()})

def scan_log(n, A, nmin, terms):
    best = None
    for al in np.linspace(0.3, 1.2, 9001):
        cols = [np.ones_like(n)] + [f(n, al) for f in terms]
        M = np.array(cols).T; y = np.log(A) + (1 + al) * np.log(n)
        s = np.linalg.lstsq(M[n >= nmin], y[n >= nmin], rcond=None)
        r = np.sum((M[n >= nmin] @ s[0] - y[n >= nmin]) ** 2)
        if best is None or r < best[0]: best = (r, al)
    return best[1]

def lokal(n, A):
    return [(int(n[i]), -np.log(A[i + 1] / A[i]) / np.log(n[i + 1] / n[i]) - 1) for i in range(len(n) - 1)]

T1 = [lambda n, a: n ** (-a)]
T2 = [lambda n, a: n ** (-a), lambda n, a: 1 / n]
T3 = [lambda n, a: n ** (-a), lambda n, a: 1 / n, lambda n, a: n ** (-2 * a)]
d8 = json.load(open(H.parent.parent / "runde_2026_10_18" / "02_theta_hoehere_raenge" / "daten" / "theta8_20_0.03_64.json"))
n8 = np.array([z[0] for z in d8["zeilen"]], float); A8 = np.array([z[1] for z in d8["zeilen"]])
print("\nRang 8 (n <= 24), lokaler Exponent - 1:", [(k, round(v, 4)) for k, v in lokal(n8, A8)][-4:])
for nmin in (3, 5, 8):
    print(f"  nmin={nmin}: alpha (T1) = {scan_log(n8, A8, nmin, T1):.4f}, (T2) = {scan_log(n8, A8, nmin, T2):.4f}, (T3) = {scan_log(n8, A8, nmin, T3):.4f}   Ziel {vor[8]:.4f}")

d10 = json.load(open(H / "daten" / "theta10_20_0.03_64.json"))
a = np.array([complex(s.replace(" ", "").replace("j", "j")) for s in d10["a"]])
x10 = complex(d10["x8"].replace(" ", "")).real
n10 = np.arange(1, 32, dtype=float); a10 = a[1:32]
A10 = np.abs(a10)
print("\nRang 10 (n <= 31), lokaler Exponent - 1:", [(k, round(v, 4)) for k, v in lokal(n10, A10)][-4:])
for nmin in (3, 5, 8):
    print(f"  nmin={nmin}: alpha (T1) = {scan_log(n10, A10, nmin, T1):.4f}, (T2) = {scan_log(n10, A10, nmin, T2):.4f}, (T3) = {scan_log(n10, A10, nmin, T3):.4f}   Ziel {vor[10]:.4f}")
# komplexes Modell
b = a10 * np.exp(2j * np.pi * n10 * x10)
for nmin in (3, 6):
    best = None
    for al in np.linspace(0.5, 1.1, 6001):
        M = np.array([np.ones_like(n10), n10 ** (-al), 1 / n10, n10 ** (-2 * al)]).T.astype(complex)
        y = b * n10 ** (1 + al); m = n10 >= nmin
        s = np.linalg.lstsq(M[m], y[m], rcond=None); r = np.sum(np.abs(M[m] @ s[0] - y[m]) ** 2)
        if best is None or r < best[0]: best = (r, al, s[0])
    print(f"  komplex nmin={nmin}: alpha = {best[1]:.4f}, c1 = {best[2][0]:.5f}, c2 = {best[2][1]:.4f}   Ziel {vor[10]:.4f}")
# Vergleich: Modell mit festem Vorhersage-alpha, Restgroesse
for al, name in ((vor[10], "Vorhersage"), (1.0, "alpha=1 (Gegenprobe)"), (0.73, "alpha=0.73 (Gegenprobe)")):
    M = np.array([np.ones_like(n10), n10 ** (-al), 1 / n10, n10 ** (-2 * al)]).T.astype(complex)
    y = b * n10 ** (1 + al); m = n10 >= 4
    s = np.linalg.lstsq(M[m], y[m], rcond=None); r = np.abs(M[m] @ s[0] - y[m]) / np.abs(y[m])
    print(f"  Rang 10, {name}: max rel. Rest (n>=4) = {r.max():.2e}")
