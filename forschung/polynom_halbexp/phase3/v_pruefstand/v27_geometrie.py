"""Strang V — geometrischer Faktor in Ks K1″(c): G_N(t) ≤ |t|·Q_k·√(N+1)·‖Ĵ_N⁻¹‖, Q_k = sup |w − w'|/|k(w) − k(w')| über 𝔻̄,
dazu c_k = (sup|k'|/2π)^{1/2}/min_γ|t| aus K1″(b). k(w) = ρw·exp(F(w)) aus der Theodorsen-Reihe (v12), F reell, |F_n| ≤ 1e-8 für n ≥ 2000.
Q_k: Schätzung über min|k'| auf ∂𝔻 (Minimumprinzip) und über Zufallspaare in 𝔻̄ (Untergrenze für das Supremum). Aufruf: python v27_geometrie.py S"""
import math
import sys

import numpy as np

S = sys.argv[1]
sys.argv = ["v12", "kubisch", "100", S]
import v12_konform_S as V12  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
F, it, d = V12.theodorsen()
F = F[:4000].real
n = np.arange(len(F))
rho = V12.rho_f


def k(w):
    return rho * w * np.exp(np.polynomial.polynomial.polyval(w, F))


def kp(w):
    Fv = np.polynomial.polynomial.polyval(w, F)
    Fp = np.polynomial.polynomial.polyval(w, F[1:] * n[1:])
    return rho * np.exp(Fv) * (1 + w * Fp)


psi = 2 * np.pi * np.arange(20000) / 20000
wb = np.exp(1j * psi)
kpb = np.abs(kp(wb))
tb = np.abs(k(wb))
rng = np.random.default_rng(1)
r1 = np.sqrt(rng.random(200000)); r2 = np.sqrt(rng.random(200000))
w1 = r1 * np.exp(2j * np.pi * rng.random(200000)); w2 = r2 * np.exp(2j * np.pi * rng.random(200000))
q = np.max(np.abs(w1 - w2) / np.abs(k(w1) - k(w2)))
qb = 1 / kpb.min()
Qk = max(q, qb)
ck = math.sqrt(kpb.max() / (2 * math.pi)) / tb.min()
print(f"S = {S}: min_∂𝔻|k'| = {kpb.min():.5f}, max|k'| = {kpb.max():.5f}; Q_k ≥ max(1/min|k'|, Zufallspaare) = max({qb:.5f}, {q:.5f}); "
      f"max_γ|t| = {tb.max():.5f}, min_γ|t| = {tb.min():.5f}; Faktor |t|·Q_k ≤ {tb.max() * Qk:.4f}; c_k = {ck:.4f}")
