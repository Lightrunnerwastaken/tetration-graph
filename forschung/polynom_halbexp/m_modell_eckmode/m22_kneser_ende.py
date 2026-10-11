"""M22 — Endstück für Kneser (Lemma M8-K): ε_L ≥ ‖1_{t≤−T1}(K̃ − A)1_{t≤−T1}‖, A = Eckmodell mit μ = √L (Symbol wie m16 'sqrtL').
Wie m17, aber die Abweichungen von G, W_in, H, W_out kommen aus mk_ecke_kneser (S2 = sup_{D̄}|F''| aus (K1)+(K2)), nicht aus dem
Schwarz-Lemma (f_K ist in L nicht analytisch). Quelle der Koeffizienten: Vorzertifikat oder Datei von Strang P.
Aufruf: python m22_kneser_ende.py "T1,..." [quelle]
"""
import sys

from flint import arb

from mk_ecke_kneser import eckschranken
from mk_kneser_arb import KneserArb

sys.stdout.reconfigure(encoding="utf-8")
T1s = [arb(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "4,5,6").split(",")]
quelle = sys.argv[2] if len(sys.argv) > 2 else "vorzertifikat"
K = KneserArb(quelle)
pi = arb.pi()
ell = K.ell
b = K.mu.arg()
amu = abs(K.mu)
print(f"Quelle: {quelle} (zertifiziert = {K.zertifiziert}); S2 ≤ {float(K.S2):.5f}; μ = √L, |μ| = {float(amu):.6f}, b = {float(b):.6f}")


def J(y_lo, y_hi, shift_sign, du, kern):
    """Wie m17: obere Schranke für ∫ e^{−Δ/2} m(Δ) dΔ (Majorante von P0/P1/|∇P0|/|∇P1| bei Höhe y ∈ [y_lo, y_hi],
    horizontale Lage shift_sign·Re ℓ ± du); ∫e^{−α/2}/(cosh α − cos y)dα = π/sin(y/2)."""
    ed = du.exp()
    vor = (-shift_sign * ell.real / 2).exp() * (du / 2).exp()
    if kern in ("P0", "G0"):
        cmax = arb(max(y_lo.cos().upper(), y_hi.cos().upper()))
        ceff = ed * cmax if bool(cmax > 0) else cmax / ed
    else:
        cmin = arb(min(y_lo.cos().lower(), y_hi.cos().lower()))
        ceff = -(cmin / ed if bool(cmin > 0) else cmin * ed)
    assert bool(ceff < 1)
    yeff = ceff.acos()
    integral = pi / (yeff / 2).sin()
    if kern in ("P0", "P1"):
        smax = arb(1) if bool(y_lo < pi / 2) and bool(y_hi > pi / 2) else arb(max(y_lo.sin().upper(), y_hi.sin().upper()))
        return vor * ed * smax / (2 * pi) * integral
    return vor * ed / (2 * pi) * integral


for T1 in T1s:
    q1 = (-T1).exp()
    e = eckschranken(K, q1)
    du0, du1 = e["gG"], e["gH"]
    yl0, yh0 = b - du0, b + du0
    yl1, yh1 = pi - b - du1, pi - b + du1
    dr = 2 * q1 / (1 - q1)
    d_in = dr * (1 / amu + e["gWin"]) + e["gWin"]
    d_out = dr * (amu + e["gWout"]) + e["gWout"]
    Jm = {}
    for ee, (yl, yh, sgn, du, dvor, wmod) in enumerate([(yl0, yh0, 1, du0, d_in, 1 / amu), (yl1, yh1, -1, du1, d_out, amu)]):
        for ep, (kP, kG) in enumerate([("P0", "G0"), ("P1", "G1")]):
            Jm[(ee, ep)] = dvor * J(yl, yh, sgn, du, kP) + wmod * du * J(yl, yh, sgn, du, kG)
    zeilen = [Jm[(0, 0)] + Jm[(0, 1)], Jm[(1, 0)] + Jm[(1, 1)]]
    spalten = [Jm[(0, 0)] + Jm[(1, 0)], Jm[(0, 1)] + Jm[(1, 1)]]
    epsL = (arb(max(z.upper() for z in zeilen)) * arb(max(z.upper() for z in spalten))).sqrt()
    print(f"T1 = {T1}: q1 = {float(q1.upper()):.3e}; |G−t−ℓ| ≤ {float(du0.upper()):.2e}, |W_in−1/μ| ≤ {float(e['gWin'].upper()):.2e}, "
          f"|H−…| ≤ {float(du1.upper()):.2e}, |W_out−μ| ≤ {float(e['gWout'].upper()):.2e}; ε_L ≤ {float(epsL.upper()):.5f}")
