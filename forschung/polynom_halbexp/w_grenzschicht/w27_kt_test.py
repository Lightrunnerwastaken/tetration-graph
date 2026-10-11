"""W27 — Welche Menge K_T trägt den Eckteil? Test: λ_A(η) (gemessen, asymptotisch, w24/w25) gegen Vorhersage ω^{K}_η(λ_D) − g_K(η) für
Zwei-Kreis-Mengen K_b = D̄ ∪ D̄(c, R) (zweiter Kreis durch L, L̄ und den reellen Scheitel b·ρ), Randdaten λ_D auf Γin (gemessen), 0 auf dem
zweiten Bogen. b = ã₀ + R̃₀ = 1.3456 ist W̄. Ebenso Vorhersage auf dem Schlitz gegen gemessen.
Aufruf: python w27_kt_test.py
"""
import json
import sys

import numpy as np

import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
NS = [500, 600, 700, 800, 900, 1000]
D = {N: json.load(open(M.W.HIER + f"/w24_N{N}.json")) for N in NS}
L, rho = M.L, M.rho
names = list(D[NS[0]]["pts"].keys())


def lnK(N, nm):
    p = D[N]["pts"][nm]
    S = D[N]["S"]
    r = max(1, sum(1 for s in S[:len(p["lnK"])] if s > 1.40))
    return 0.5 * np.log(sum(np.exp(2 * v) for v in p["lnK"][:r]))


Ns = np.array(NS, float)
A = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T
lam = {}
for nm in names:
    co, *_ = np.linalg.lstsq(A, np.array([lnK(N, nm) for N in NS]), rcond=None)
    lam[nm] = co[2]
for nm in names:
    if nm != "L":
        lam[nm] -= lam["L"]
lam["L"] = 0.0
z = {nm: complex(*D[NS[0]]["pts"][nm]["z"]) for nm in names}


def make_map(b):
    """τ-Abbildung für K_b (unskaliert): zweiter Kreis durch L, L̄, b·ρ."""
    B = b * rho
    th1 = np.angle((-rho - L) / (-rho - np.conj(L)))
    th2 = np.angle((B - L) / (B - np.conj(L)))

    def inside(lo, hi, t):
        return ((t - lo) % (2 * np.pi)) < ((hi - lo) % (2 * np.pi))
    if inside(th2, th1, 0.0):
        lo, wid = th2, (th1 - th2) % (2 * np.pi)
    else:
        lo, wid = th1, (th2 - th1) % (2 * np.pi)
    tau0 = np.exp(1j * ((0.0 - lo) % (2 * np.pi)) * np.pi / wid)

    def tau(zeta):
        m = (zeta - L) / (zeta - np.conj(L))
        ang = (np.angle(m) - lo) % (2 * np.pi)
        return (abs(m) ** (np.pi / wid)) * np.exp(1j * min(ang, wid) * np.pi / wid), (ang <= wid)
    return tau, tau0, wid


# D-Bogen-Profil (gemessen) als Funktion des Winkels τ°; Spiegelung symmetrisch
dpts = sorted([(float(nm[1:]), lam[nm]) for nm in names if nm[0] == "D"])
dt = np.array([0.0] + [a for a, _ in dpts])
dl = np.array([0.0] + [b for _, b in dpts])
argL = np.angle(L)


def lamD(phi):            # phi = Winkel auf |x| = ρ (rad), Γin: argL ≤ |phi| ≤ π
    t = np.degrees(np.abs(np.angle(np.exp(1j * phi))) - argL)
    return np.interp(np.clip(t, 0, 103.38), dt, dl)


for b in (1.3456, 1.30, 1.27, 1.24, 1.20, 1.15):
    tau, t0, wid = make_map(b)
    # Randdaten in s: Γin-Punkte über Winkel, zweiter Bogen = 0
    phis = np.linspace(argL, 2 * np.pi - argL, 4001)[1:-1]
    sv = np.array([tau(rho * np.exp(1j * p))[0].real for p in phis])
    lv = lamD(phis)
    o = np.argsort(sv)
    sv, lv = sv[o], lv[o]

    def omega(zeta):
        tz, ok = tau(zeta)
        if not ok:
            return np.nan, np.nan
        x, y = tz.real, tz.imag
        # Poisson-Integral über die Γin-Halbachse (Stützstellen sv), Rest 0
        ker = y / ((sv - x) ** 2 + y ** 2) / np.pi
        om = np.trapezoid(lv * ker, sv)
        g = np.log(abs((tz - np.conj(t0)) / (tz - t0)))
        return om, g
    rows = []
    for nm in names:
        if nm[0] in "AS":
            om, g = omega(z[nm])
            rows.append((nm, lam[nm], om - g if not np.isnan(om) else np.nan))
    err = [abs(a - b_) for _, a, b_ in rows if not np.isnan(b_)]
    print(f"b={b:.4f} (Aussenwinkel {np.degrees(wid):.1f}°): mittlere |λ_gem − (ω−g)| = {np.mean(err):.4f}")
    print("   " + " ".join(f"{nm}:{a:+.4f}/{b_:+.4f}" for nm, a, b_ in rows))
