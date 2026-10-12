"""Strang V — Satz V2 mit strengen Green-Konstanten (Auftrag Leitung 02:05, Punkt 2): Bilanz wie v31, aber
  g_γ(L(1+u)) ≥ g_{K*}(L(1+u))   (Lune K* = D₁ ∪ D₂ ⊃ K_γ, v39, streng)
  dist(L(1+u), K_γ) ≥ dist(L(1+u), K*) = min_i(|η − c_i| − r_i)   (streng, explizit)
  Kontur Γ*_d = {g_{K*} = d}: dort g_γ ≥ d (streng); Länge aus v39 (obere Schranke); Schlitz bis u₀(d) ≤ 0.2
  Lip(g_γ) ≤ 0.8606 (v40, innere Scheiben)          M_γ ≤ 2.40 (v32)
  M_d = max|f_K| auf Γ*_d: f(x) = exp(f(log x)) (log x ∈ D̄), f(x) = log f(eˣ) (eˣ ∈ D̄) mit Kugeln (mk_kneser_arb); Punkte nahe dem Schlitz-
       Schnitt, wo keine der beiden Formeln greift: lokale Schranke |f(x)| ≤ |L| + |μ||x − L| + ½C₂|x − L|² mit C₂ = 8 (v32: sup|F''| am Segment ≤ 7.99/ρ … konservativ).
Aufruf: python v41_N0_streng.py S d
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from scipy.integrate import quad

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
from mk_kneser_arb import KneserArb  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
S = sys.argv[1]
d = float(sys.argv[2])
import os
LU = json.loads((HIER / f"v39_lune_S{float(S)}_t{os.environ.get('V_TAU1', '12')}_p{os.environ.get('V_P', '1.3')}.json").read_text(encoding="utf-8"))
assert LU["linse_ok"] and LU["fout_ok"]
KN = KneserArb(str(HIER / "kneser_halbexp_b_zert_V.json"))
ctx.prec = 120
rho, th = float(KN.rho.mid()), float(KN.th.mid())
L = complex(rho * math.cos(th), rho * math.sin(th))
c1, c2 = LU["c1"], LU["c2"]
r1, r2 = abs(c1 - L), abs(c2 - L)
lip, Mg = float(os.environ.get("V_LIP", "0.8606")), 2.40                   # Keil c = 0.25: V_LIP = 3.2196 (v40)
CJ, SIG, K0, c0 = float(os.environ.get("V_CJ", "0.95")), float(os.environ.get("V_SIG", "5.447")), 6.0, 0.641   # streng: V_CJ = 1.8595, V_SIG = 5.4469 bis u ≤ 0.35 (v44)
UJ = float(os.environ.get("V_UJ", "0.2"))
minfp, max_m_fp, M2 = 0.34558, 2.11804, 0.81835          # streng (v43)
tab = np.array(LU["schlitz"])
uu, gl = tab[:, 0], tab[:, 1]
niv = LU["niveau"][str(d)]
u0, laenge, pts = niv["u0"], niv["laenge"], [complex(a, b) for a, b in niv["punkte"]]
assert u0 <= UJ, f"u₀ = {u0} > {UJ}"


def g_u(u):
    """Treppen-Unterschranke: g_γ(L(1+u)) ≥ g_{K*}(L(1+u_k)) für u ≥ u_k (g_{K*} längs des Schlitzes monoton steigend, in v39 geprüft)."""
    k = np.searchsorted(uu, u, side="right") - 1
    return 0.0 if k < 0 else float(gl[k])


def dist_u(u):
    eta = L * (1 + u)
    return max(min(abs(eta - c1) - r1, abs(eta - c2) - r2), 0.0)


# M_d streng (Kugeln) bzw. lokal
Lb = acb(KN.rho * KN.th.cos(), KN.rho * KN.th.sin())
werte, lokal = [], 0
zweifach = [0]
# Kugeln über die Bögen zwischen benachbarten Stützpunkten (Radius = halbe Sehne + 0.01 für die Bogenabweichung)
kugeln = []
maxsehne = max(abs(pts[(i + 1) % len(pts)] - pts[i]) for i in range(len(pts)))
assert maxsehne < 0.1, f"Sehne {maxsehne} zu gross"
for i in range(len(pts)):
    x0, x1 = pts[i], pts[(i + 1) % len(pts)]
    m = (x0 + x1) / 2
    r = abs(x1 - x0) / 2 + 0.7 * abs(x1 - x0) ** 2 / 8 + 1e-6          # Bogenabweichung ≤ κ·Sehne²/8, κ ≤ 0.7 (Γ*_d explizit)
    kugeln.append(acb(arb(m.real, r), arb(m.imag, r)))
for X in kugeln:
    lx = X.log()
    try:
        if bool(abs(lx) < KN.rho * arb("0.999")):
            v = (KN.F_kugel(lx / KN.rho) * KN.rho).exp()
        elif bool((X.exp() / KN.rho).abs_upper() < 0.999):
            v = (KN.F_kugel(X.exp() / KN.rho) * KN.rho).log()
        elif bool(abs(lx.log()) < KN.rho * arb("0.999")):
            # nahe dem Schlitz-Schnitt: zweimal log zieht zu L hin (|log'(L)| = 1/ρ, Drehung −θ), f = exp∘exp∘f∘log∘log
            v = (KN.F_kugel(lx.log() / KN.rho) * KN.rho).exp().exp()
            zweifach[0] += 1
        else:
            raise AssertionError
        werte.append(float(abs(v).upper()))
    except AssertionError:
        dx = abs(X - Lb)
        v = abs(Lb) + abs(KN.mu) * dx + 4 * dx ** 2
        werte.append(float(v.upper()))
        lokal += 1
Md = max(werte)
print(f"[S = {S}, d = {d}] α = {LU['alpha']:.5f}, u₀ = {u0:.4f}, |Γ*_d| ≤ {laenge:.3f}, M_d ≤ {Md:.3f} ({zweifach[0]} Bogenkugeln über log∘log, {lokal} von {len(kugeln)} lokal/offen), "
      f"Lip ≤ {lip}, M_γ ≤ {Mg}", flush=True)
rmin = rho + max(niv["min_abs"] - rho, 1e-3)


def A_theta(N):
    return (5 / rho) / minfp + max_m_fp / (rho * c0) + M2 / (c0 * N)


def terme(N):
    sq = math.sqrt(N + 1)
    ePre = math.exp((N + 2) * lip / N)

    def EB(u):
        dist = dist_u(u)
        if dist >= 1.0 / N:
            return sq * N * (1 + Mg) * ePre * math.exp(-(N + 2) * g_u(u))
        return sq * (1 + Mg) / max(dist, 1e-300)

    def EC(u):
        dist = dist_u(u)
        ea = rho * (1 + u)
        if dist >= 1.0 / N:
            return (N + 1) / rho * ePre * math.exp(-(N + 2) * g_u(u)) + (rho / ea) ** (N + 1) / (ea - rho)
        return 1.0 / (ea - rho)

    pre = 2 * rho / (2 * math.pi)
    uk = 1 / (N * 0.7 * rho)
    A_s = pre * quad(lambda u: CJ * u ** SIG * EB(u), 0, u0, limit=500, points=[uk])[0]
    C_s = pre * quad(lambda u: CJ * u ** SIG * EC(u), 0, u0, limit=500, points=[uk])[0]
    A_c = Md * laenge / (2 * math.pi) * sq * N * (1 + Mg) * ePre * math.exp(-(N + 2) * d)
    C_c = Md * laenge / (2 * math.pi) * ((N + 1) / rho * ePre * math.exp(-(N + 2) * d) + (rho / rmin) ** (N + 1) / (rmin - rho))
    return 2 * (A_s + A_c), C_s + C_c


bmax = {}
for N in [300, 318, 350, 400, 500, 700, 1000, 2000, 5000]:
    A, C = terme(N)
    b = (0.5 / (A_theta(N) * N * K0) - C) / A
    bmax[N] = b
    print(f" N = {N:4d}: A_N = {A:.3e}, C_N = {C:.3e}, B_max(N) = {b:.4g}", flush=True)
for N0 in (318, 400):
    print(f"⇒ B_max(N₀ = {N0}) = {min(v for k, v in bmax.items() if k >= N0):.4g}  (streng bis auf: C_J, A_θ-Randwerte, Lemmata X5/X7/X8′/X4a, (Fam))", flush=True)
