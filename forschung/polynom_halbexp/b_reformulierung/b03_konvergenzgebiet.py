"""b03: Konvergenzgebiet der p_N in der Ebene. Karte von log10|p_N - p_194| (N = 40, 80, 120) und log10|p_N - f_K|,
f_K ausserhalb D über f(z) = exp(f(log z)) bzw. f(z) = log f(e^z). Ausgabe: Textkarte + PNG."""
import json, sys
import numpy as np
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
L = complex(-mp.lambertw(-1)); L = L.conjugate() if L.imag < 0 else L

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

xs = np.linspace(-2.6, 2.6, 53)
ys = np.linspace(0.0, 2.6, 27)
NS = (40, 80, 120, 160)
D = {N: np.zeros((len(ys), len(xs))) for N in NS}
for i, y in enumerate(ys):
    for j, x in enumerate(xs):
        z = mp.mpc(x, y)
        ref = horner(P[194], z)
        for N in NS:
            D[N][i, j] = float(mp.log10(abs(horner(P[N], z) - ref) + mp.mpf(10) ** -38))
np.savez("b03_karte.npz", xs=xs, ys=ys, **{f"N{N}": D[N] for N in NS})

# Textkarte: Konvergenzrate zwischen N=80 und 160 (Stellen pro 80 Grad, negativ = Konvergenz)
print("Rate log10|p_160 - p_194| - log10|p_80 - p_194| (negativ = konvergiert). Zeilen y von oben.")
for i in reversed(range(0, len(ys), 2)):
    row = ""
    for j in range(0, len(xs), 2):
        r = D[160][i, j] - D[80][i, j]
        ch = "#" if r < -6 else ("+" if r < -2 else ("." if r < 1 else "X"))
        if abs(complex(xs[j], ys[i])) < abs(L) and abs(complex(xs[j], ys[i]) - L) > 0.08:
            ch = ch.lower() if ch != "." else ","
        row += ch
    print(f"y={ys[i]:4.1f} {row}")
print("Legende: # >6 Stellen besser, + 2..6, ./, -2..1 (kaum), X divergiert; Kleinbuchstaben/Komma = innerhalb D")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
for a, N in zip(ax, (80, 160)):
    im = a.imshow(D[N], origin="lower", extent=(xs[0], xs[-1], ys[0], ys[-1]), vmin=-30, vmax=5, cmap="viridis")
    a.add_patch(plt.Circle((0, 0), abs(L), fill=False, color="w"))
    a.plot([L.real], [L.imag], "r*")
    a.set_title(f"log10|p_{N} - p_194|")
    fig.colorbar(im, ax=a)
fig.savefig("b03_karte.png", dpi=110, bbox_inches="tight")
print("Bild: b03_karte.png")
