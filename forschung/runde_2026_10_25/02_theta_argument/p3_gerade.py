"""Gebiet 02, Teil E: Φ über eine ganze Periode als Phase von h' auf der Geraden Im w = π.

Exakt (Kettenregel, Hauptzweige): für s ∈ (0,1) ist sexp(s−2) = −X < 0 und
   Φ(s) = arg(1+θ'(s)) = arg h'(sexp(s−2) + i0) = P(Y) − π,   P(Y) := arg h'(Y + iπ),   Y = ln(−sexp(s−2)).
s → 1⁻  ⇔  Y → −∞ ;   s → 0⁺  ⇔  Y → +∞.  Asymptotisch arg a_n ≈ π + Φ(σ_n), σ_n → 0⁺, also wandert arg a_n wie P(ln ln n).
Ausgabe: stetiger Zweig von P längs der ganzen Geraden (Y = ±e^u), Spannweite, Abgleich mit Teil B.
Aufruf: python p3_gerade.py
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from kern import Fix

HIER = Path(__file__).resolve().parent
mp.mp.dps = 30
fx = Fix(mp.e, K=40)
lam = fx.lam


def hp(w):
    w = mp.mpc(w)
    acc = mp.mpc(1)
    n = 0
    while abs(w - fx.L) > mp.mpf("1e-3") or n < 5:
        acc /= w
        w = mp.log(w)
        n += 1
    u = w - fx.L
    tol = mp.mpf(10) ** (-mp.mp.dps)
    while abs(u) ** (fx.K + 1) > tol * abs(u):
        acc /= (fx.L + u)
        u = mp.log(fx.L + u) - fx.L
    S = fx.S
    su = mp.polyval(S[::-1], u)
    dsu = mp.polyval([k * S[k] for k in range(len(S) - 1, 0, -1)], u)
    return acc * dsu / su / mp.log(lam)


def P(Y):
    return float(mp.arg(hp(mp.mpc(Y, mp.pi))))


# Abgleich mit Teil B (s aus Newton, x = sexp(s)):  Y = ln(−ln ln x)
print("Abgleich: Φ(s) aus Teil B gegen P(Y) − π")
for x, phiB in [("1.2", 0.5777476), ("1.5", 0.1462485), ("2.0", -0.3532388), ("2.3", -0.6131248)]:
    Y = mp.log(-mp.log(mp.log(mp.mpf(x))))
    v = P(Y) - np.pi
    v = (v + np.pi) % (2 * np.pi) - np.pi
    print(f"  x={x}: Y={float(Y):+.6f}  P(Y)−π = {v:+.7f}   Teil B: {phiB:+.7f}")

# stetiger Zweig längs der Geraden, Y = ∓e^u; adaptive Verfeinerung bis jeder Phasenschritt < 0.3; Normierung Φ(1/2-Bereich) nahe 0
U = 19.5


def spur(u0, u1, sgn, n0):
    grid = list(np.linspace(u0, u1, n0))
    out = [(grid[0], P(sgn * mp.exp(grid[0])))]
    stack = [(u, None) for u in reversed(grid[1:])]
    while stack:
        u, val = stack.pop()
        if val is None:
            val = P(sgn * mp.exp(u))
        d = (val - out[-1][1] + np.pi) % (2 * np.pi) - np.pi
        if abs(d) > 0.3 and abs(u - out[-1][0]) > 1e-12:
            stack.append((u, val))
            stack.append((0.5 * (u + out[-1][0]), None))
            continue
        out.append((u, val))
    return out


neg = spur(U, -4.0, -1, 1500)
pos = spur(-4.0, U, +1, 1500)
seq = [(-u, v) for u, v in neg] + [(0.0, P(0))] + [(u, v) for u, v in pos]
uu = np.array([q[0] for q in seq])
ph = np.unwrap(np.array([q[1] for q in seq]))
k0 = len(neg)
Phi = ph - np.pi
Phi -= 2 * np.pi * np.round(Phi[k0] / (2 * np.pi))
print("Punkte:", len(seq))
print(f"\nΦ längs der Periode (stetiger Zweig, Φ(Y=0) = {Phi[k0]:+.6f}):")
for tag, sel in [("Y<0 (s→1⁻)", slice(0, k0)), ("Y>0 (s→0⁺)", slice(k0 + 1, None))]:
    seg = Phi[sel]
    print(f"  {tag}: min {seg.min():+.5f}  max {seg.max():+.5f}")
print(f"  gesamt: min {Phi.min():+.5f}  max {Phi.max():+.5f}  Spannweite {Phi.max()-Phi.min():.5f}")
# ausgewählte Stellen (ln|Y| = Turmwerte usw.)
probe = [0, 1, np.e, 3.5, 5, 10, np.e ** np.e, 16, 19]
print("\n ln|Y| | Φ(Y<0)     | Φ(Y>0)")
rows = []
for p in probe:
    i_neg = int(np.argmin(abs(uu[:k0] + p))) if p > 0 else k0
    i_pos = k0 + 1 + int(np.argmin(abs(uu[k0 + 1:] - p))) if p > 0 else k0
    rows.append([float(p), float(Phi[i_neg]), float(Phi[i_pos])])
    print(f"{p:6.3f} | {Phi[i_neg]:+.6f} | {Phi[i_pos]:+.6f}")
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / "p3_e.json").write_text(json.dumps({"probe": rows, "min": float(Phi.min()), "max": float(Phi.max()),
                                                     "u_step": 0.004, "U": U}, indent=1), encoding="utf-8")
np.save(HIER / "daten" / "p3_phi.npy", np.vstack([uu, Phi]))
