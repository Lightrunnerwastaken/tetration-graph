"""W37 — Green-Funktion von Ĉ∖K_T, K_T = W̄ ∪ ℬ (Linse aus w35/w36, x₀ = 0.7), per Ladungssimulation; Vergleich mit g_W̄ auf dem Schlitz
und auf Γ_d (W̄): ε-Werte für Satz W1 (Bernstein–Walsh auf K_T). Unskaliert.
"""
import sys

import numpy as np

import w_kern as W
import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
exec(open(W.HIER + "/w36_satzW1_pruefung.py", encoding="utf-8").read().split("ZT = [")[0].split("import w24_rand_genau as M")[1])
rho, L = W.RHO, W.L
argL = np.angle(L)
a0, R0 = M.a0, M.R0
psiL = np.angle(L - a0)
# Rand von K_T (obere Hälfte, dann gespiegelt): Linsenrand Z⁺ von L bis −ρ(1+s), A-Bogen von L bis a₀+R₀
u = np.linspace(0, 1, 500)
taus_z = (180 - np.degrees(argL)) * (np.sinh(4 * u) / np.sinh(4))       # dicht bei L
Zp = np.array([rho * np.exp(s_lens(argL + np.deg2rad(t)) + 1j * (argL + np.deg2rad(t))) for t in taus_z[1:]])
psA = psiL * (1 - np.sinh(4 * u) / np.sinh(4))                            # von psiL (L) bis 0 (Scheitel), dicht bei L
Ap = a0 + R0 * np.exp(1j * psA[1:])
oben = np.concatenate([Zp[::-1], [L], Ap])                                # von −ρ über L zum A-Scheitel
rand = np.concatenate([oben, np.conj(oben[::-1])[1:-1]])                  # geschlossene Kurve
# Ladungen nach innen versetzt
t_ = np.gradient(np.concatenate([rand[-1:], rand, rand[:1]]))[1:-1]
nr = 1j * t_ / np.abs(t_)
mid = 0.5 + 0j
if np.mean(np.abs(rand + 0.01 * nr - mid)) > np.mean(np.abs(rand - 0.01 * nr - mid)):
    nr = -nr
ds = np.abs(t_)
src = rand + nr * np.minimum(4 * ds, 0.5 * np.abs(rand - L) * np.abs(rand - np.conj(L)) + 1e-9)
Mat = np.log(np.abs(rand[:, None] - src[None, :]))
Mat = np.vstack([np.hstack([Mat, np.ones((len(rand), 1))]), np.hstack([np.ones((1, len(src))), np.zeros((1, 1))])])
rhs = np.concatenate([np.zeros(len(rand)), [1.0]])
sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=1e-14)
q, C = sol[:-1], sol[-1]
print(f"Kollokationsrest {np.abs((Mat @ sol - rhs)[:-1]).max():.1e}")


def gKT(t):
    return float(np.sum(q * np.log(np.abs(t - src))) + C)


print("Schlitz L(1+u): g_W̄, g_{K_T}, Verhältnis")
for uu in (1e-3, 3e-3, 0.01, 0.03, 0.06, 0.1, 0.151):
    z = L * (1 + uu)
    gw, gk = M.green(z), gKT(z)
    print(f"  u = {uu}: {gw:.5f}  {gk:.5f}  {gk/gw:.3f}")
for dd in (0.05, 0.1):
    rat = [gKT(z) / dd for z in M.gamma_d_punkte(dd, 37)]
    print(f"Γ_{dd}: min g_KT/d = {min(rat):.3f}  (⇒ ε_Γ ≥ {min(rat):.3f})")
print(f"Probe: g_KT(3) = {gKT(3+0j):.4f}, g_W̄(3) = {M.green(3+0j):.4f}")
