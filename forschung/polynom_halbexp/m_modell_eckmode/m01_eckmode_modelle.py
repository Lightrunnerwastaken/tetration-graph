"""M01 — Welche Modelle haben die Eckmode? σ₁..σ₃ der ρ-skalierten J_N⁻¹ (ℓ²_ρ).
Aufruf: python m01_eckmode_modelle.py modell N1,N2,... [param]
  modell: kubisch (param = μ-Ausdruck, Standard sqrt(L)), quad (param = α), affin (μ), kquad (komplex: L + μ(x−L) + b(x−L)², param "mu,b")
"""
import sys
import time

import numpy as np

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
modell = sys.argv[1]
NS = [int(x) for x in sys.argv[2].split(",")]
par = sys.argv[3] if len(sys.argv) > 3 else None
L = K.L_KOMPLEX
B = int(sys.argv[4]) if len(sys.argv) > 4 else None

if modell == "kubisch":
    mu = eval(par or "np.sqrt(L)", {"np": np, "L": L})
    coef = K.kubisch(mu)
elif modell == "quad":
    coef = K.quadratisch(float(par or 0.2283))
elif modell == "affin":
    mu = eval(par or "np.sqrt(L)", {"np": np, "L": L})
    coef = [L * (1 - mu), mu]
elif modell == "kquad":
    mu, b = [eval(t, {"np": np, "L": L}) for t in par.split(";")]
    # f = L + μ(x−L) + b(x−L)²
    coef = [L - mu * L + b * L * L, mu - 2 * b * L, b]
else:
    raise SystemExit("unbekanntes Modell")

print(f"Modell {modell}, Koeffizienten {coef}", flush=True)
print(f"f(L) − L = {abs(K.polywert(coef, L) - L):.1e}, f'(L) = {K.polyabl(coef, L):.5f} (|·| = {abs(K.polyabl(coef, L)):.5f})", flush=True)
phi = np.linspace(0, 2 * np.pi, 7201)
x = K.RHO * np.exp(1j * phi)
fx = K.polywert(coef, x)
inn = np.abs(fx) < K.RHO
print(f"Anteil Γin: {inn.mean():.3f}; max|f| auf |x|=ρ: {np.abs(fx).max():.4f}; f(0) = {coef[0]}", flush=True)
t0 = time.time()
for N in NS:
    S, U, Vt, cols = K.ecksingulaer(coef, N, B=B)
    Fu = np.abs(K.auf_kreis(U[:, 0]))
    M = len(Fu)
    print(f"N={N} ({time.time()-t0:.0f}s): σ = {', '.join(f'{s:.5f}' for s in S)} | u1 max bei {360*int(np.argmax(Fu))/M:.2f}°",
          flush=True)
