"""D06 — Das Symbol von G = lim J_N^{-1}: Vergleich der Toeplitz-Koeffizienten im Inneren mit
       σ(φ) = f'(x) e^{−x} · 1[|f(x)| < |L|],   x = |L| e^{iφ}   (f ≈ p_N),
       g_d = (1/2π) ∫ σ(φ) e^{−idφ} dφ   (Eintrag (j+d, j) der umskalierten Matrix |L|^{k−j} (J^{-1})_{kj}).
Liest d05_N<N>.npz (aus d05_struktur.py).
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
L = 1.3745570107437075
Ji = np.load(f"d05_N{N}.npz")["Ji"]
n = N + 1
k = np.arange(n)
G = Ji * (L ** k)[:, None] / (L ** k)[None, :]
a = [float(x) for x in json.load(open("../gemeinsam/polynome_bis_194.json"))[str(N)]]
pp = np.polynomial.Polynomial(a)
dp = pp.deriv()
M = 1 << 16
phi = 2 * np.pi * np.arange(M) / M
x = L * np.exp(1j * phi)
inn = np.abs(pp(x)) < L
sig = dp(x) * np.exp(-x) * inn
gd_theorie = np.fft.fft(sig) / M          # gd_theorie[d] = (1/M) Σ σ e^{−i d φ}
D = 30
fehler = []
print(" d   Matrix (Mittel j=50..90)     Theorie σ·1_in            1/m ohne Abschneiden")
sig_voll = dp(x) * np.exp(-x)
gd_voll = np.fft.fft(sig_voll) / M
for d in list(range(-6, 7)) + [10, -10, 20, -20, 30, -30]:
    vals = np.array([G[j + d, j] for j in range(50, 91)])
    m = vals.mean()
    th = gd_theorie[d % M]
    fehler.append(abs(m - th))
    print(f"{d:3d}  {m:+.5f} (±{vals.std():.1e})   {th.real:+.5f}{th.imag:+.5f}i   {gd_voll[d % M].real:+.5f}{gd_voll[d % M].imag:+.5f}i")
print("max |Matrix − Theorie| über diese d:", f"{max(fehler):.2e}")
phiL = np.angle(complex(0.31813150520476413, 1.3372357014306895))
print("Sprung bei arg L =", f"{np.degrees(phiL):.3f}°;  Rand von {{|f(x)|<|L|}} auf dem Kreis bei",
      f"{np.degrees(phi[np.argmax(inn[:M//2])]):.3f}°", "; 1/m(L) = 1/sqrt(L) =", 1 / np.sqrt(complex(0.31813150520476413, 1.3372357014306895)))
