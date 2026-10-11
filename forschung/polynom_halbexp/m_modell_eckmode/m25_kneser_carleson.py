"""M25 — Carleson-Konstante C_h und m_in für KNESER (wie M20; Kanten aus mk_kneser_arb, Ecken aus mk_ecke_kneser).
C_h² ≤ (r_max/d0)·[(1+ε) + (1+1/ε)(y2−y1)²κ], Höhe y = Im H ∈ [y1, y2], d0 ≤ Re H', r = |ξ'(t+iπ)|/|ξ'(H)|.
H'(t) = [1/(v−L) − 1/(v−L̄)]·y'(t)/F'(v), y = ξ(t+iπ), v = h(y), y'(t) = (L − L̄)M/(M−1)², M = −e^{iθ}e^t.
m_in = min_Γin |m|, |m(z)| = e^{ρ Re z}/|F'(z)|. Symmetrie t ↦ −t (reelle Koeffizienten).
Aufruf: python m25_kneser_carleson.py [Ta] [Zellweite] [c] [quelle]
"""
import sys

import numpy as np
from flint import acb, arb

from mk_ecke_kneser import eckschranken
from mk_kneser import h_newton
from mk_kneser_arb import KneserArb

sys.stdout.reconfigure(encoding="utf-8")
Ta = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(13)
dz = arb(sys.argv[2]) if len(sys.argv) > 2 else arb("0.02")
c_akk = arb(sys.argv[3]) if len(sys.argv) > 3 else arb("0.8")
quelle = sys.argv[4] if len(sys.argv) > 4 else "vorzertifikat"
KN = KneserArb(quelle)
pi = arb.pi()
th = KN.th
L = KN.eith
Lb = L.conjugate()
emith = acb(th.cos(), -th.sin())


def xi(M):
    return Lb + (Lb - L) / (M - 1)


def winkel(x):
    x0 = acb(x.real.mid(), x.imag.mid())
    r = x0 / abs(x0)
    return (x * r.conjugate()).arg() + r.arg()


def absxi_abl(zeta):
    M = L * zeta.exp()
    return 2 * th.sin() * zeta.real.exp() / abs(M - 1) ** 2


def daten(t):
    M = -L * t.exp()
    y = xi(M)
    ph = winkel(y)
    y0 = complex(float(y.real.mid()), float(y.imag.mid()))
    v0 = complex(h_newton(np.array([y0]))[0])
    v = KN.h(acb(ph.cos(), ph.sin()), v0)
    assert v is not None
    H = ((v - L) / (v - Lb) * emith).log()
    yp = (L - Lb) * M / (M - 1) ** 2
    Hp = (1 / (v - L) - 1 / (v - Lb)) * yp / KN.Fp_kugel(v)
    r = absxi_abl(acb(t, pi)) / absxi_abl(H)
    return H, Hp.real, r


N = int(float((Ta / dz).upper()) + 0.5)
y1 = y2 = d0 = None
rmax = arb(0)
for i in range(N):
    m = -Ta + dz * (arb(i) + arb("0.5"))
    tb = arb(m.mid(), (dz / 2).upper() * (1 + 1e-12))
    H, hp, r = daten(tb)
    y1 = H.imag.lower() if y1 is None else arb(min(y1, H.imag.lower()))
    y2 = H.imag.upper() if y2 is None else arb(max(y2, H.imag.upper()))
    d0 = hp.lower() if d0 is None else arb(min(d0, hp.lower()))
    rmax = arb(max(rmax, r.upper()))
print(f"Quelle {quelle} (zertifiziert = {KN.zertifiziert}); Raster [−{Ta}, 0]: Im H ∈ [{float(y1):.5f}, {float(y2):.5f}], "
      f"Re H' ≥ {float(d0):.5f}, r ≤ {float(rmax):.5f}")
E = eckschranken(KN, (-Ta).exp())
b = E["b"]
y1e, y2e = pi - b - E["gH"], pi - b + E["gH"]
print(f"Ecke t < −{Ta}: Im H ∈ [{float(y1e.lower()):.5f}, {float(y2e.upper()):.5f}], Re H' ≥ {float(E['reHp'].lower()):.5f}, "
      f"r ≤ {float(E['rmax'].upper()):.5f}")
Y1 = arb(min(y1, y1e.lower()))
Y2 = arb(max(y2, y2e.upper()))
D0 = arb(min(d0, E["reHp"].lower()))
RM = arb(max(rmax, E["rmax"].upper()))
kap = arb(max((1 / (arb(1).exp() * (pi - Y2))) ** 2, (1 / (arb(1).exp() * Y1)) ** 2).upper())
best = None
for k in range(1, 200):
    e_ = arb(k) / 1000
    val = (RM / D0) * ((1 + e_) + (1 + 1 / e_) * (Y2 - Y1) ** 2 * kap)
    if best is None or val.upper() < best[0].upper():
        best = (val, e_)
Ch = best[0].sqrt()
# m_in über Γin = {θ ≤ |φ| ≤ π}: |m| = e^{ρ cos φ}/|F'(e^{iφ})|
nphi = 4000
m_in = None
for k in range(nphi):
    a_ = th + (pi - th) * arb(k) / nphi
    b_ = th + (pi - th) * arb(k + 1) / nphi
    ph = arb(((a_ + b_) / 2).mid(), ((b_ - a_) / 2).upper() * (1 + 1e-12))
    v = ((KN.rho * ph.cos()).exp() / abs(KN.Fp_kreis(ph))).lower()
    m_in = v if m_in is None else arb(min(m_in, v))
CE = (1 / m_in ** 2 + Ch ** 2).sqrt()
print(f"κ ≤ {float(kap):.4f}: C_h ≤ {float(Ch.upper()):.5f}; min_Γin|m| ≥ {float(m_in):.6f}; C_E ≤ {float(CE.upper()):.5f}")
print(f"(a″)_∞ Kneser: σ_min(J_∞|P_N) ≥ c/C_E ≥ {float((c_akk / CE).lower()):.4f} (mit c = {c_akk} aus der IMS-Zusammenstellung)")

# --- Verfeinerung (Kanten getrennt): ‖𝒯x‖² ≤ ‖Z0‖²/m_in² + (r/d0)∫|Z(s+iy(s))|²ds und für y1 = min Höhe, θy = y1/π:
# ‖Z(·+iy1)‖² = ∫|Ẑ|²e^{−2σy1} ≤ α‖Z0‖² + β‖Z1‖² mit α = (1−θy)s^θy, β = θy s^{θy−1} (Young, x = e^{−2σ}, beliebiges s > 0);
# Ableitungsterm ≤ κ'‖Z‖², κ' = (1+1/ε)(y2−y1)²κ. ⇒ C_E'² = max(1/m_in² + R[(1+ε)α + κ'], R[(1+ε)β + κ']), R = r_max/d0.
def ce_fein(m_in, RM, D0, Y1, Y2, kap):
    R = RM / D0
    th_y = Y1 / pi
    beste = None
    for ke in range(1, 120):
        e_ = arb(ke) / 1000
        kp = (1 + 1 / e_) * (Y2 - Y1) ** 2 * kap
        for js in range(1, 200):
            s_ = arb(js) / 200
            al = (1 - th_y) * s_ ** th_y
            be = th_y * s_ ** (th_y - 1)
            A = 1 / m_in ** 2 + R * ((1 + e_) * al + kp)
            Bv = R * ((1 + e_) * be + kp)
            v = arb(max(A.upper(), Bv.upper()))
            if beste is None or v < beste[0]:
                beste = (v, e_, s_)
    return beste


bf = ce_fein(m_in, RM, D0, Y1, Y2, kap)
CEf = bf[0].sqrt()
print(f"Verfeinert (Kanten getrennt, ε = {float(bf[1]):.3f}, s = {float(bf[2]):.3f}): C_E' ≤ {float(CEf.upper()):.5f} ⇒ σ_min(J_∞|P_N) ≥ "
      f"{float((c_akk / CEf).lower()):.4f}")
