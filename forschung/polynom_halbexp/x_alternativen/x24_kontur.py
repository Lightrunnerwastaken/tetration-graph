"""X24 — Konturkonstanten für Satz P-I im W̄-Rahmen: Niveaulinie Γ_d = {g_W̄ = d} (Green-Funktion der Zwei-Kreis-Menge explizit über Möbius +
Potenz), Länge |Γ_d|, Abstand δ_d = dist(Γ_d, W̄), max|Γ_d| (Lage gegenüber den nächsten Singularitäten log log a₀ = −0.362 ± πi, |x| = 3.16)
und M_d = max_{Γ_d}|f_K| (f_K ausserhalb D über f(w) = exp(f(log w)); Punkte, wo die Rekursion scheitert, werden gezählt und ausgelassen —
das sind die Schlitzkreuzungen). Ausgabe der Konturkonstante M_d|Γ_d|/(2πδ_d) für d = 0.1, 0.2, 0.3, 0.5.
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
mp.mp.dps = 30
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(L)
a0 = mp.mpf(A[0])
R0 = abs(L - a0)
Am = [mp.mpf(s) for s in A[:1500]]


def fK(z, tiefe=0):
    if abs(z) < 0.97 * rho:
        s = mp.mpc(0)
        for c in reversed(Am):
            s = s * z + c
        return s
    if tiefe > 7:
        raise RuntimeError
    return mp.exp(fK(mp.log(z), tiefe + 1))


Lc, Lb = complex(L), complex(mp.conj(L))
rho_f, a0_f, R0_f = float(rho), float(a0), float(R0)


def inW(z):
    return abs(z) <= rho_f + 1e-12 or abs(z - a0_f) <= R0_f + 1e-12


# Möbius T(z) = (z − L)/(z − L̄): Aussengebiet → Sektor; Öffnung und Richtung aus Randpunkten bestimmen
T = lambda z: (z - Lc) / (z - Lb)
# Randbögen von ∂W̄ ausserhalb des jeweils anderen Kreises
ph = np.linspace(0, 2 * np.pi, 200001)
c1 = rho_f * np.exp(1j * ph)
c1 = c1[np.abs(c1 - a0_f) > R0_f]
c2 = a0_f + R0_f * np.exp(1j * ph)
c2 = c2[np.abs(c2) > rho_f]
a1 = np.angle(T(c1[len(c1) // 2]))
a2 = np.angle(T(c2[len(c2) // 2]))
tinf = np.angle(T(1e9))   # Bild von ∞ = 1 ⇒ Winkel 0
# Sektor zwischen a1 und a2, der tinf enthält
lo, hi = sorted([a1, a2])
if not (lo < tinf < hi):
    lo, hi = hi, lo + 2 * np.pi
open_ = hi - lo
alpha = open_ / np.pi
print(f"Sektoröffnung {np.degrees(open_):.3f}° (α = {alpha:.6f}; Soll 0.882964)")


def g(z):
    t = T(z)
    ang = np.angle(t)
    if ang < lo:
        ang += 2 * np.pi
    s = (np.abs(t) ** (1 / alpha)) * np.exp(1j * (ang - lo) / alpha)   # obere Halbebene
    p = 1.0 ** (1 / alpha) * np.exp(1j * (tinf - lo if tinf >= lo else tinf - lo + 2 * np.pi) / alpha)
    return np.log(np.abs((s - np.conj(p)) / (s - p)))


for d in (0.1, 0.2, 0.3, 0.5):
    pts = []
    for th in np.linspace(0, 2 * np.pi, 721)[:-1]:
        # Strahl vom Zentrum 0.5 aus; Bisektion auf g = d
        c0 = 0.6 + 0j
        rlo, rhi = 0.0, 8.0
        e = np.exp(1j * th)
        while not inW(c0 + rlo * e) or False:
            break
        # finde ersten Punkt ausserhalb W̄
        r = 0.0
        while inW(c0 + r * e):
            r += 0.01
        rlo = r
        while g(c0 + rhi * e) < d:
            rhi *= 1.5
        for _ in range(60):
            mid = 0.5 * (rlo + rhi)
            if g(c0 + mid * e) < d:
                rlo = mid
            else:
                rhi = mid
        pts.append(c0 + rlo * e)
    P = np.array(pts)
    laenge = np.sum(np.abs(np.diff(np.append(P, P[0]))))
    # Abstand zu W̄
    bd = np.concatenate([c1[::50], c2[::50]])
    dist = min(np.abs(P[:, None] - bd[None, :]).min(axis=1))
    Mmax, fehl = 0.0, 0
    for z in P[::4]:
        try:
            Mmax = max(Mmax, float(abs(fK(mp.mpc(z.real, z.imag)))))
        except (RuntimeError, ValueError, ZeroDivisionError):
            fehl += 1
    print(f"d = {d}: |Γ_d| = {laenge:.3f}, δ_d = {dist:.4f}, max|z| auf Γ_d = {np.abs(P).max():.3f}, min|z − (−0.362 ± πi)| = "
          f"{min(np.abs(P - (-0.3624 + 3.1416j)).min(), np.abs(P - (-0.3624 - 3.1416j)).min()):.3f}, M_d = {Mmax:.3f} ({fehl} Punkte ohne Wert); "
          f"Konturkonstante M_d|Γ_d|/(2πδ_d) = {Mmax * laenge / (2 * np.pi * dist):.2f}", flush=True)
