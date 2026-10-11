"""W06 — Phasenraumbilder (Gabor-Transformation im Index k, Winkel φ) der Eckmoden, kubisches Modell.
|Σ_k a_k g(k − k0) e^{−ikφ}| für Fensterzentren k0/N ∈ [0,1], φ ∈ [0°, 360°); log10, normiert.
Aufruf: python w06_phasenraum_bild.py N   → w06_phasenraum_N{N}.png
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(W.HIER + f"/w01_N{N}.npz")
S, U, Vt = d["S"], d["U"], d["Vt"]
n = N + 1
k = np.arange(n)
argL = np.degrees(np.angle(W.L))
M = 2048
x0s = np.linspace(0.0, 1.0, 201)
wbr = 0.6 * np.sqrt(N)


def portrait(a):
    img = np.zeros((len(x0s), M))
    for i, x0 in enumerate(x0s):
        g = np.exp(-0.5 * ((k - x0 * N) / wbr) ** 2)
        seq = np.zeros(M, dtype=complex)
        seq[:n] = a * g
        # Σ_k a_k g e^{−ikφ}: fft liefert Σ seq_k e^{−2πi k m/M}
        img[i] = np.abs(np.fft.fft(seq))
    return np.log10(img + 1e-300)


fig, axs = plt.subplots(2, 2, figsize=(14, 10))
for ax, (name, a) in zip(axs.flat, [("v1 (Eingabe/adjungiert)", Vt[0]), ("v2", Vt[1]), ("u1 (Ausgabe)", U[:, 0]), ("u2", U[:, 1])]):
    P = portrait(a)
    P -= P.max()
    phi = np.arange(M) * 360 / M
    im = ax.imshow(P, origin="lower", aspect="auto", extent=[0, 360, 0, 1], vmin=-12, vmax=0, cmap="viridis")
    ax.axvline(argL, color="r", lw=0.5)
    ax.axvline(360 - argL, color="r", lw=0.5)
    ax.set_title(f"N={N}: {name}")
    ax.set_xlabel("φ (Grad)")
    ax.set_ylabel("k/N")
    plt.colorbar(im, ax=ax)
plt.tight_layout()
plt.savefig(W.HIER + f"/w06_phasenraum_N{N}.png", dpi=80)
print("ok")
