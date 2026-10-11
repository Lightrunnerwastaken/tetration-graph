"""W32 — Selbstkonsistentes Leitordnungsproblem (n-te Wurzel) für die adjungierte Eckmode — Versuch, c und λ herzuleiten.

Befunde (w29–w31): H := Φ_y (innen, Ω₁) bzw. Φ_T (aussen, Ĉ∖D̄) ist EINE harmonische Funktion auf 𝒟' = Ĉ∖(Bucht̄), Bucht = Gebiet zwischen
Γout und der Kurve B; H ~ −log|t| bei 0 und ∞; H = 0 auf B; und im Halbmond gilt T ≈ 𝒴, also H(t) = H(h(t)) (h = Hauptzweig von F̃⁻¹).
Ansatz H = −log|t| + v, v beschränkt harmonisch auf 𝒟' (Ladungen in der Bucht, Σq = 0):
   v = log|t| auf B,   v(t) − v(h(t)) = −log|h(t)| für t ∈ Γout (|t| = 1).
Freier Rand: B ← h(Z_T), Z_T = {H = 0} im Halbmond. Iteration ab B₀ = h(A-Bogen). Ausgabe: c = H(−1), Profil H(e^{iφ}) auf Γin, Buchtrand bei 0°.
Rechnung umskaliert (t = x/ρ). Aufruf: python w32_selbstkonsistent.py [Iterationen] [mu_abs mu_arg]
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
ITER = int(sys.argv[1]) if len(sys.argv) > 1 else 8
mu = float(sys.argv[2]) * np.exp(1j * float(sys.argv[3])) if len(sys.argv) > 3 else np.sqrt(W.L)
coef, Lf = W.kubik_koeff(mu=mu)
rho = W.RHO
Ft = np.polynomial.Polynomial([coef[j] * rho ** (j - 1) for j in range(4)])
Lt = Lf / rho
phiL = np.angle(Lt)
a0 = Ft(0).real
R0 = abs(Lt - a0)


def h(t):
    """Hauptzweig: Wurzel von F̃(z) = t mit kleinstem |z − t·0.6|… stetig: nimm die Wurzel mit kleinstem |z|."""
    r = (Ft - t).roots()
    return r[np.argmin(np.abs(r))]


hv = np.vectorize(h)

# Startkurve B₀ = h(A-Bogen): Polardarstellung r_B(φ)
psis = np.linspace(np.angle(Lt - a0), -np.angle(Lt - a0), 2001)
A = a0 + R0 * np.exp(1j * psis)
B0 = hv(A)
phB = np.angle(B0)
o = np.argsort(phB)
phi_grid = np.linspace(-phiL, phiL, 401)
rB = np.interp(phi_grid, phB[o], np.abs(B0)[o])
rB[0] = rB[-1] = 1.0


def loese(rB):
    nB = len(phi_grid)
    Bpts = rB * np.exp(1j * phi_grid)
    # Kollokation auf B (ohne Endpunkte) und auf Γout
    Bc = Bpts[2:-2:2]
    Gc = np.exp(1j * np.linspace(-phiL, phiL, 241)[3:-3])
    hG = hv(Gc)
    # Ladungen in der Bucht: zwischen B und Γout, mehrere Radien, mit Abstand zu beiden Rändern
    srcs = []
    for frac in (0.15, 0.5, 0.85):
        rr = rB + frac * (1 - rB)
        srcs.append(rr[3:-3:2] * np.exp(1j * phi_grid[3:-3:2]))
    src = np.concatenate(srcs)
    def Lmat(pts):
        return np.log(np.abs(pts[:, None] - src[None, :]))
    MB = np.hstack([Lmat(Bc), np.ones((len(Bc), 1))])
    MG = np.hstack([Lmat(Gc) - Lmat(hG), np.zeros((len(Gc), 1))])
    Msum = np.hstack([np.ones((1, len(src))), np.zeros((1, 1))])
    Mat = np.vstack([MB, MG, 1e3 * Msum])
    rhs = np.concatenate([np.log(np.abs(Bc)), -np.log(np.abs(hG)), [0.0]])
    sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=1e-13)
    res = Mat @ sol - rhs
    q, C = sol[:-1], sol[-1]

    def H(t):
        t = np.atleast_1d(t)
        return -np.log(np.abs(t)) + np.log(np.abs(t[:, None] - src[None, :])) @ q + C
    return H, np.abs(res[:len(Bc)]).max(), np.abs(res[len(Bc):len(Bc) + len(Gc)]).max()


for it in range(ITER):
    H, rb, rg = loese(rB)
    c = H(np.array([-1.0 + 0j]))[0]
    # Z_T im Halbmond: auf Strahlen ψ (|ψ| < phiL) den Radius r > 1 mit H = 0 (Bisektion)
    Z = []
    for ps in np.linspace(-phiL * 0.995, phiL * 0.995, 161):
        lo, hi = 1.0 + 1e-6, 3.0
        flo = H(np.array([lo * np.exp(1j * ps)]))[0]
        if flo <= 0:
            continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if H(np.array([mid * np.exp(1j * ps)]))[0] > 0:
                lo = mid
            else:
                hi = mid
        Z.append(0.5 * (lo + hi) * np.exp(1j * ps))
    Z = np.array(Z)
    Bn = hv(Z)
    ph = np.angle(Bn)
    o = np.argsort(ph)
    rBn = np.interp(phi_grid, ph[o], np.abs(Bn)[o], left=1.0, right=1.0)
    rBn[0] = rBn[-1] = 1.0
    print(f"It {it}: Rest B {rb:.1e}, Γout {rg:.1e} | c = H(−1) = {c:+.5f} | Bucht(0°) = {np.interp(0, phi_grid, rB):.4f} → "
          f"{np.interp(0, phi_grid, rBn):.4f} | Z_T(0°) = {abs(Z[len(Z)//2]):.4f}", flush=True)
    rB = 0.5 * rB + 0.5 * rBn
H, rb, rg = loese(rB)
print(f"Ende: c = {H(np.array([-1.0+0j]))[0]:+.5f}")
print("Profil H(e^{i(argL+τ)}): " + " ".join(f"{t}:{H(np.array([np.exp(1j*(phiL+np.deg2rad(t)))]))[0]:+.4f}" for t in (1, 2, 4, 8, 15, 30, 60, 103.38)))
hh = 1e-5
print("Stegfrequenz −∂_rH: " + " ".join(f"{p}°:{-(H(np.array([(1+hh)*np.exp(1j*np.deg2rad(p))]))[0]-H(np.array([(1-hh)*np.exp(1j*np.deg2rad(p))]))[0])/(2*hh):.4f}" for p in (80, 90, 120, 180)))
print("H innen (r = 0.3/0.5/0.7, φ = 0/90/180): " + " ".join(f"{H(np.array([r*np.exp(1j*np.deg2rad(p))]))[0]:.4f}" for r in (0.3, 0.5, 0.7) for p in (0, 90, 180)))
