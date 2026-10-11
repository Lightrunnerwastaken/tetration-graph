"""W35 — Konstanten für Satz W1 (kubisch, aus w01): (H-V1) Kreisprofil der Eckmoden Y_i = σ_iV_i (alle σ_i > 1.40), (H-I) Indexschwanz,
Linse ℬ = {|t| ≥ ρ: P_ext[Λ](t) > x₀ log(|t|/ρ)}, deren Green-Breite M_ℬ = max_ℬ g_W̄ und ε_Γ = 1 − M_ℬ/d, Schlitzbedingung.
Λ(φ) = asymptotisches Profil (Fit a + κ ln N + λN über N) + Fehler, auf Γin; 0 auf Γout.
Aufruf: python w35_satzW1_konstanten.py 500,600,700,800,900,1000 [x0]
"""
import sys

import numpy as np

import w_kern as W
import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
X0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7
rho, L = W.RHO, W.L
argL = np.angle(L)
taus = np.concatenate([np.geomspace(0.05, 10, 60), np.linspace(10.5, 180 - np.degrees(argL), 120)])   # Grad ab arg L, bis 180°
lnY = {}
idx = {}
for N in NS:
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, Vt = d["S"], d["Vt"]
    r = max(1, sum(1 for s in S[:Vt.shape[0]] if s > 1.40))
    z = rho * np.exp(1j * (argL + np.deg2rad(taus)))
    Y = np.sqrt(sum(S[i] ** 2 * np.abs(W.V_auswerten(Vt[i], z)) ** 2 for i in range(r)))
    lnY[N] = np.log(Y)
    yk = np.sqrt(sum((S[i] * Vt[i]) ** 2 for i in range(r)))
    idx[N] = {x: yk[: int(x * N)].max() for x in (0.5, 0.6, 0.7, 0.75, 0.8)}
Ns = np.array(NS, float)
A = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T
Ymat = np.array([lnY[N] for N in NS])
co = np.linalg.lstsq(A, Ymat, rcond=None)[0]
co1 = np.linalg.lstsq(A[:-1], Ymat[:-1], rcond=None)[0]
co2 = np.linalg.lstsq(A[1:], Ymat[1:], rcond=None)[0]
lam = co[2]
err = np.maximum(np.abs(co1[2] - lam), np.abs(co2[2] - lam))
# Ecke: λ(τ→0) als Referenz (Eckwert polynomial): λ_rel = λ − λ(τ_min)
lam_rel = lam - lam[0]
Lam_env = np.maximum(lam_rel + err, 0.0)
print("(H-V1) asymptotisches Kreisprofil Λ (λ + Fehler, relativ zur Ecke), Auswahl:")
for t in (0.1, 0.5, 1, 2, 4, 8, 15, 30, 60, 103.38):
    j = np.argmin(abs(taus - t))
    print(f"   τ={taus[j]:7.2f}°: λ = {lam_rel[j]:+.5f} ± {err[j]:.5f}  κ = {co[1][j]:+.2f}")
print(f"   λ_max = {Lam_env.max():.5f}")
# Exponent nahe L aus Λ_env
sel = (taus > 0.9) & (taus < 4.5) & (Lam_env > 0)
pp = np.polyfit(np.log(taus[sel]), np.log(Lam_env[sel]), 1)[0]
print(f"   Exponent p (Fit 1°–4.5°) = {pp:.2f}")
print("(H-I) max_{k < xN} |y_k| (alle Eckmoden):")
for x in (0.5, 0.6, 0.7, 0.75, 0.8):
    vals = [idx[N][x] for N in NS]
    rate = np.polyfit(Ns, np.log(vals), 1)[0]
    kap = np.polyfit(np.log(Ns), np.log(vals), 1)[0]
    print(f"   x = {x}: " + " ".join(f"{v:.3g}" for v in vals) + f"   (Rate {rate:+.5f}/Grad bzw. N^{kap:.2f})")


# Λ als Funktion des Winkels φ (rad, ganzer Kreis, symmetrisch)
def Lam_of_phi(phi):
    t = np.degrees(np.abs(np.angle(np.exp(1j * phi)))) - np.degrees(argL)
    return np.where(t > 0, np.interp(t, taus, Lam_env, left=0.0, right=Lam_env[-1]), 0.0)


PH = np.linspace(-np.pi, np.pi, 8001)[:-1]
LV = Lam_of_phi(PH)


def P_ext(t):
    """Poisson-Integral des Äusseren der Kreisscheibe |x| = ρ für Randdaten Λ."""
    r = abs(t) / rho
    ph = np.angle(t)
    ker = (r ** 2 - 1) / (r ** 2 - 2 * r * np.cos(ph - PH) + 1)
    return np.mean(LV * ker)


print(f"\nLinse ℬ für x₀ = {X0}: Breite s(φ) (|t| = ρe^s) mit P_ext[Λ] = x₀ s; g_W̄ auf dem Linsenrand:")
Mb = 0.0
rows = []
for t in (0.5, 1, 2, 4, 8, 15, 30, 60, 103.38):
    ph = argL + np.deg2rad(t)
    lo, hi = 0.0, 0.3
    if P_ext(rho * np.exp(1e-7 + 1j * ph)) <= X0 * 1e-7:
        s = 0.0
    else:
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if P_ext(rho * np.exp(mid + 1j * ph)) > X0 * mid:
                lo = mid
            else:
                hi = mid
        s = hi
    g = M.green(rho * np.exp(s + 1j * ph))
    Mb = max(Mb, g)
    rows.append((t, s, g))
print("   " + " ".join(f"{t}°: s={s:.5f}, g={g:.5f}" for t, s, g in rows))
# feiner Scan für M_ℬ
for t in np.linspace(0.2, 180 - np.degrees(argL), 200):
    ph = argL + np.deg2rad(t)
    lo, hi = 0.0, 0.3
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if P_ext(rho * np.exp(mid + 1j * ph)) > X0 * mid:
            lo = mid
        else:
            hi = mid
    Mb = max(Mb, M.green(rho * np.exp(hi + 1j * ph)))
print(f"   M_ℬ = max g_W̄ auf ℬ = {Mb:.5f}  ⇒  ε_Γ = 1 − M_ℬ/d = {1 - Mb/0.05:.3f} (d = 0.05), {1 - Mb/0.1:.3f} (d = 0.1)")
# Bedingung auf Γout-Seite und Schlitz: C₁ := sup (1/s) P_ext[Λ](ρe^{s+iφ}) für kleine s, φ ∈ Γout ∪ {arg L}
C1 = max(P_ext(rho * np.exp(1e-4 + 1j * ph)) / 1e-4 for ph in np.linspace(-argL, argL, 201))
print(f"   C₁ (Γout-Seite, s → 0) = {C1:.4f}  (Bedingung C₁ < x₀ = {X0})")
for u in (1e-3, 1e-2, 0.05, 0.15):
    t = L * (1 + u)
    print(f"   Schlitz u = {u}: P_ext[Λ] = {P_ext(t):.5f} gegen x₀ ln|t|/ρ = {X0*np.log(abs(t)/rho):.5f}  ⇒  {'ausserhalb ℬ' if P_ext(t) <= X0*np.log(abs(t)/rho) else 'IN ℬ'}")
np.savez(W.HIER + "/w35_profil.npz", taus=taus, Lam_env=Lam_env, lam=lam_rel, err=err)
