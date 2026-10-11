"""W24 — Randdaten für das Poisson–Jensen-Kriterium, GENAU (arb, ctx.cap gesetzt), kubisches Modell.
Für Punkte η auf ∂W̄ (D-Bogen, A-Bogen; obere Hälfte), in L, auf dem Schlitz L(1+u) und auf Γ_d (d = 0.05, 0.1):
  y_η := Π_N k_η − T_N k_η = J⁻¹T_N J∞k_η − T_N k_η  (Cs kopffreie Pol-Antwort),
  R_N(η) = Σ_j |y_j| ρ^j (ℓ¹_ρ),  K_i(η) = ⟨u_i, D y_η⟩ (Eckanteile, u_i aus w01, umskaliert), K_N = (Σ_{i≤r}|K_i|²)^{1/2}.
Ausgabe: w24_N{N}.json mit ln R_N, ln|K_1|, ln K_N(r) für r = 1..6 (soweit vorhanden), Punktlisten.
Aufruf: python w24_rand_genau.py N1,N2,...
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho = W.RHO
L = W.L
a0 = coef[0]
R0 = abs(L - a0)
psiL = np.angle(L - a0)
argL = np.angle(L)


# ---------- Green-Abbildung von Ĉ∖W̄ (unskaliert): τ(ζ) in die obere Halbebene, L ↦ 0, L̄ ↦ ∞ ----------
def tau_map(zeta):
    m = (zeta - L) / (zeta - np.conj(L))
    th1 = np.angle((-rho - L) / (-rho - np.conj(L)))
    th2 = np.angle((a0 + R0 - L) / (a0 + R0 - np.conj(L)))

    def inside(lo, hi, t):
        return ((t - lo) % (2 * np.pi)) < ((hi - lo) % (2 * np.pi))
    if inside(th2, th1, 0.0):
        lo, wid = th2, (th1 - th2) % (2 * np.pi)
    else:
        lo, wid = th1, (th2 - th1) % (2 * np.pi)
    ang = (np.angle(m) - lo) % (2 * np.pi)
    tau = (abs(m) ** (np.pi / wid)) * np.exp(1j * min(ang, wid) * np.pi / wid)
    tau0 = np.exp(1j * ((0.0 - lo) % (2 * np.pi)) * np.pi / wid)
    return tau, tau0, wid


def green(zeta):
    t, t0, _ = tau_map(zeta)
    return float(np.log(abs((t - np.conj(t0)) / (t - t0))))


def schlitz_punkt(u):
    return L * (1 + u)


def gamma_d_punkte(dd, n_pts=19):
    """Punkte auf {g = d} (obere Hälfte): Strahlen vom Punkt 0.25 aus, Bisektion im Radius."""
    out = []
    for th in np.linspace(0, np.pi, n_pts):
        lo_r, hi_r = 0.3, 6.0
        e = np.exp(1j * th)
        # innerster Punkt ausserhalb W̄ auf dem Strahl: g<0 innen (tau_map liefert dort min(ang, wid) → g≈0); robust: Bisektion auf g − d
        for _ in range(60):
            mid = 0.5 * (lo_r + hi_r)
            z = 0.25 + mid * e
            inW = (abs(z) <= rho) or (abs(z - a0) <= R0)
            gv = -1.0 if inW else green(z)
            if gv < dd:
                lo_r = mid
            else:
                hi_r = mid
        out.append(0.25 + 0.5 * (lo_r + hi_r) * e)
    return out


TAUS = [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 11, 15, 20, 30, 45, 60, 80, 103.38]
ADIST = [0.002, 0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.6, 1.9]
USL = [1e-3, 3e-3, 0.01, 0.03, 0.06, 0.1, 0.151]


def a_bogen(dist):
    dpsi = 2 * np.arcsin(min(dist / (2 * R0), 1.0))
    return a0 + R0 * np.exp(1j * (psiL - dpsi))


def punkte():
    P = [("L", L)]
    P += [(f"D{t}", rho * np.exp(1j * (argL + np.deg2rad(t)))) for t in TAUS]
    P += [(f"A{d}", a_bogen(d)) for d in ADIST]
    P += [(f"S{u}", schlitz_punkt(u)) for u in USL]
    for dd in (0.05, 0.1):
        P += [(f"G{dd}_{i}", z) for i, z in enumerate(gamma_d_punkte(dd))]
    return P


if __name__ == "__main__":
    PTS = punkte()
    for N in [int(x) for x in sys.argv[1].split(",")]:
        n = N + 1
        ctx.prec = int(0.9 * N) + 400
        ctx.cap = n + 10
        t0 = time.time()
        d = np.load(W.HIER + f"/w01_N{N}.npz")
        U, S = d["U"], d["S"]
        a = W.kubik_a(N, coef)
        P = W.potenzen(a, N)
        J = W.J_mat(a, N, P)
        x0 = acb_poly([acb(c) for c in a])
        mser = acb_poly([0])
        for j in range(1, n):
            mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
        mser = mser.truncate(n)
        E = arb_mat(n, 2 * len(PTS))
        KZ = []
        for i, (nm, z) in enumerate(PTS):
            zz = acb(z.real, z.imag)
            s1 = list(acb_series((acb_poly([zz]) - x0).coeffs(), prec=n).inv().coeffs())
            s1 += [acb(0)] * (n - len(s1))
            gl = [zz ** (-(k + 1)) for k in range(n)]
            s2 = list((mser * acb_poly(gl)).truncate(n).coeffs())
            s2 += [acb(0)] * (n - len(s2))
            for k in range(n):
                w = s1[k] + s2[k]
                E[k, 2 * i] = w.real
                E[k, 2 * i + 1] = w.imag
            KZ.append(gl)
        X = J.solve(E, algorithm="approx")
        rr = arb(rho)
        Dk = [rr ** k for k in range(n)]
        uarb = [[arb(float(U[k, c])) for k in range(n)] for c in range(U.shape[1])]
        out = {"S": [float(s) for s in S[:8]], "pts": {}}
        for i, (nm, z) in enumerate(PTS):
            y = [acb(X[k, 2 * i], X[k, 2 * i + 1]) - KZ[i][k] for k in range(n)]
            yD = [y[k] * Dk[k] for k in range(n)]
            Rn = sum((abs(v) for v in yD), arb(0))
            Ks = []
            for c in range(U.shape[1]):
                s = acb(0)
                for k in range(n):
                    s += uarb[c][k] * yD[k]
                Ks.append(float(abs(s).mid().log()) if abs(s).mid() > 0 else -1e9)
            out["pts"][nm] = {"z": [z.real, z.imag], "lnR": float(Rn.mid().log()), "lnK": Ks}
        json.dump(out, open(W.HIER + f"/w24_N{N}.json", "w"))
        print(f"N={N} ({time.time()-t0:.0f}s): lnR(L) = {out['pts']['L']['lnR']:.3f}, lnK1(L) = {out['pts']['L']['lnK'][0]:.3f}, "
              f"Punkte {len(PTS)}", flush=True)
