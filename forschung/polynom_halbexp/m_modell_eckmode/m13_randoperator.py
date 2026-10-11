"""M13 — Der Randoperator 𝒯 im radialen Modell (nur Werte in D̄, keine Abschneidekante):
  (𝒯x)(φ) = x(e^{iφ}) + x(f(e^{iφ}))/m(φ)       für φ ∈ Γin   (f(e^{iφ}) = s(φ)e^{iφ} ∈ D)
  (𝒯x)(φ) = x(e^{iφ}) + f'(e^{iφ})·x(h(e^{iφ}))  für φ ∈ Γout  (h(e^{iφ}) ∈ D)
Kern von 𝒯 auf H²(D) = Kern von T auf E²(Ω). Hier: σ_min(𝒯|P_N) in L²(∂D) (Gewicht 1), float, Gitterquadratur.
Aufruf: python m13_randoperator.py "N1,N2,..." [M Gitterpunkte]
"""
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)


def f(z):
    return eps + c0 * z + eps * z * z


def fp(z):
    return c0 + 2 * eps * z


def h(w):
    # Zweig mit h(a0) = 0: εz² + c0 z + ε − w = 0, z = (−c0 + sqrt(c0² − 4ε(ε − w)))/(2ε)
    return (-c0 + np.sqrt(c0 * c0 - 4 * eps * (eps - w) + 0j)) / (2 * eps)


# Gitter (Mittelpunktsregel, je Bogen getrennt, Knicke bei ±θ sind Gitterränder)
def bogen(a, b, n):
    t = (np.arange(n) + 0.5) / n
    return a + (b - a) * t, (b - a) / n


nin = int(M * (1 - th / np.pi))
nout = M - nin
phin, win = bogen(th, 2 * np.pi - th, nin)
phout, wout = bogen(-th, th, nout)
zin, zout = np.exp(1j * phin), np.exp(1j * phout)
gin = f(zin)                       # = s e^{iφ}
min_ = fp(f(zin))                  # m = f'∘f
gout = h(zout)
fpo = fp(zout)
print(f"Kontrollen: max|f(e^iφ)| auf Γin = {np.abs(gin).max():.4f}, max|h(e^iφ)| auf Γout = {np.abs(gout).max():.4f}, "
      f"|f(h(z)) − z| max = {np.abs(f(gout) - zout).max():.1e}", flush=True)
for N in NS:
    t0 = time.time()
    k = np.arange(N + 1)
    Ain = zin[:, None] ** k[None, :] + (gin[:, None] ** k[None, :]) / min_[:, None]
    Aout = zout[:, None] ** k[None, :] + fpo[:, None] * gout[:, None] ** k[None, :]
    A = np.vstack([Ain * np.sqrt(win / (2 * np.pi)), Aout * np.sqrt(wout / (2 * np.pi))])
    S = np.linalg.svd(A, compute_uv=False)
    print(f"N={N}: σ_min(𝒯|P_N) = {S[-1]:.6f}, nächste {S[-2]:.6f}, {S[-3]:.6f}; σ_max = {S[0]:.4f} ({time.time()-t0:.0f}s)", flush=True)
