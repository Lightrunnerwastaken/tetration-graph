"""Strang C, Runde 4, c28: Konstanten der Eck-Bedingung (E) — Rate gegen Gebiet.

Geometrie (c19/c24): E = Ĉ∖W̄, τ: E → obere Halbebene, Ecke L ↦ 0, L̄ ↦ ∞, D-Bogen ↦ (−∞,0), A-Bogen ↦ (0,∞), ∞ ↦ τ0.
Für ein Ratenprofil λ auf ∂W (Funktion des Randparameters s = τ(η) ∈ ℝ) ist ω_ζ(λ) = (1/π)∫ λ(s) Im τ/|s−τ|² ds.
(E_∂) verlangt Θ(λ) := sup_{ζ∈Schlitz∪Γ'} ω_ζ(λ)/g(ζ) < 1 (Γ' = Niveaulinie {g = d}: dort einfach ω_ζ(λ) < d).
Hier: (1) Umrechnung s ↔ |η − L| ↔ Winkel; (2) Stufenprofil λ = λ0 auf dem D-Bogen ausserhalb der beiden Eckumgebungen
|s| < s1 bzw. |s| > |τ0|²/s1 (A-Bogen: λ ≤ 0, Strang D: Eckanteil fällt ausserhalb |ζ| = ρ) → kritisches λ*(s1);
(3) Grenzwert am Schlitz: lim_{u→0} ω/g = (1/κ0)∫ λ(s) s^{−2} ds, κ0 = 2π Im τ0/|τ0|²;
(4) Probe des Grenzwerts gegen die exakte Formel.
"""
import sys

sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, log, arg, conj, atan, findroot, quad, inf

mp.dps = 30
L = mpc("0.31813150520476413", "1.3372357014306895")
Lq = conj(L)
RHO = abs(L)
a0 = mpf("0.49856328794111")
R0 = abs(L - a0)


def M(z):
    return (z - L) / (z - Lq)


th_a, th_b = arg(M(-RHO)), arg(M(a0 + R0))
lo, hi = sorted([th_a, th_b])
if not (lo < 0 < hi):
    lo, hi = hi, lo + 2 * pi
ALPHA = (hi - lo) / pi


def tau(z):
    w = M(z) * exp(-1j * lo)
    a = arg(w)
    if a < 0:
        a += 2 * pi
    a = min(max(a, mpf(0)), ALPHA * pi)
    return abs(w) ** (1 / ALPHA) * exp(1j * a / ALPHA)


TAU0 = tau(mpf(10) ** 30)
KAPPA0 = 2 * pi * TAU0.imag / abs(TAU0) ** 2


def g(z):
    t = tau(z)
    return log(abs((t - conj(TAU0)) / (t - TAU0)))


def om_int(t, a, b):
    x, y = t.real, t.imag
    fa = atan((a - x) / y) if a != -inf else -pi / 2
    fb = atan((b - x) / y) if b != inf else pi / 2
    return (fb - fa) / pi


aL = arg(L)
seite = -1 if tau(-RHO * (1 + mpf(10) ** -20)).real < 0 else 1
print(f"α = {mp.nstr(ALPHA, 7)}, τ0 = {mp.nstr(TAU0, 7)}, |τ0| = {mp.nstr(abs(TAU0), 7)}, κ0 = 2π Im τ0/|τ0|² = {mp.nstr(KAPPA0, 6)},"
      f" D-Bogen auf der Seite {'τ<0' if seite < 0 else 'τ>0'}")

# (1) Umrechnung: Punkt auf dem D-Bogen im Abstand r von L  ↔  s
print("\n(1) D-Bogen: Abstand |η−L|, Winkel arg η − arg L, Randparameter s = τ(η):")
for r in (0.002, 0.005, 0.0107, 0.02, 0.05, 0.1, 0.2):
    phi = 2 * mp.asin(r / (2 * RHO))
    eta = RHO * exp(1j * (aL + phi)) * (1 + mpf(10) ** -20)
    print(f"   |η−L| = {r}: Winkel {mp.nstr(phi * 180 / pi, 4)}°, s = {mp.nstr(tau(eta).real, 5)}")


def r_von_s(s1):
    """Abstand zu L auf dem D-Bogen mit |τ| = s1."""
    f = lambda ph: abs(tau(RHO * exp(1j * (aL + ph)) * (1 + mpf(10) ** -20)).real) - s1
    ph = findroot(f, mpf(s1) ** ALPHA * 2)
    return 2 * RHO * mp.sin(ph / 2)


# Niveaulinie {g = d}: Punkte über τ (g = d ist im τ-Bild ein Kreis um τ0-Achse): τ mit |(τ−τ̄0)/(τ−τ0)| = e^d
def niveau_tau(d, n=720):
    # Apollonius-Kreis: Mittelpunkt c, Radius r
    k = exp(d)
    c = (k ** 2 * TAU0 - conj(TAU0)) / (k ** 2 - 1)
    r = k * abs(TAU0 - conj(TAU0)) / (k ** 2 - 1)
    return [c + r * exp(2j * pi * j / n) for j in range(n)]


def Theta_stufe(s1, d, nslit=60):
    """sup über Schlitz (u ≤ u0(d)) von ω(B)/g und über {g=d} von ω(B)/d, B = D-Bogen ohne Eckumgebungen."""
    a_, b_ = (-abs(TAU0) ** 2 / s1, -s1) if seite < 0 else (s1, abs(TAU0) ** 2 / s1)
    # Schlitz: u von 1e-6 bis u0 mit g(L(1+u0)) = d
    u0 = findroot(lambda u: g(L * (1 + u)) - d, mpf("0.1"))
    best = mpf(0)
    for j in range(nslit + 1):
        u = mpf(10) ** (-6 + (log(u0, 10) + 6) * j / nslit)
        z = L * (1 + u)
        best = max(best, om_int(tau(z), a_, b_) / g(z))
    best_n = mpf(0)
    for t in niveau_tau(d):
        best_n = max(best_n, om_int(t, a_, b_) / d)
    return best, best_n, u0


print("\n(2) Stufenprofil: λ = λ0 auf dem D-Bogen ausser |s| < s1 (bei L) und |s| > |τ0|²/s1 (bei L̄); A-Bogen λ ≤ 0.")
print("    kritische Rate λ*(s1) = 1/max(Θ_Schlitz, Θ_Niveau); Grenzwert-Formel κ0·s1/(1 − s1²/|τ0|²)")
for s1 in (0.001, 0.0019, 0.005, 0.0109, 0.03, 0.0663, 0.143):
    r1 = r_von_s(mpf(s1))
    zeile = []
    for d in (0.05, 0.1):
        ts, tn, u0 = Theta_stufe(mpf(s1), mpf(d))
        zeile.append(f"d={d} (u0={mp.nstr(u0, 3)}): Θ_S={mp.nstr(ts, 4)}, Θ_Γ={mp.nstr(tn, 4)}, λ*={mp.nstr(1 / max(ts, tn), 4)}")
    print(f"   s1 = {s1} (|η−L| = {mp.nstr(r1, 3)}): " + " | ".join(zeile)
          + f" | Grenzwert κ0 s1/(1−s1²) = {mp.nstr(KAPPA0 * s1 / (1 - s1 ** 2 / abs(TAU0) ** 2), 4)}")

# (3) Probe des Grenzwerts für ein glattes Profil λ(s) = λ0 |s|^p/(s2^p + |s|^p) auf dem D-Bogen (p = 2)
print("\n(3) Glattes Profil λ(s) = λ0·s²/(s2² + s²) auf dem D-Bogen (bis s = −1, dann symmetrisch s ↦ |τ0|²/s):")
for s2 in (0.002, 0.01):
    lam0 = mpf("0.012")
    prof = lambda s: lam0 * s ** 2 / (s2 ** 2 + s ** 2)
    # Grenzwert-Integral (1/κ0)∫ λ/s² ds über s ∈ (−∞, 0) mit Symmetrie: 2·∫_{−1}^{0} (bis auf ∫ über |s|>1 ~ s1)
    I = quad(lambda s: prof(s) / s ** 2, [0, s2, 1]) + quad(lambda s: prof(1 / s) / s ** 2, [1, 1 / s2, inf])
    grenz = I / KAPPA0
    u = mpf("1e-5")
    z = L * (1 + u)
    t = tau(z)
    om = quad(lambda s: prof(-s) * t.imag / abs(-s - t) ** 2, [0, s2, 1]) / pi \
        + quad(lambda s: prof(-1 / s) * t.imag / abs(-s - t) ** 2, [1, 1 / s2, inf]) / pi
    print(f"   s2 = {s2}, λ0 = 0.012: (1/κ0)∫λ/s² = {mp.nstr(grenz, 4)}, exakt ω/g bei u = 1e-5: {mp.nstr(om / g(z), 4)}")

# (4) Richtung des Schlitzes im τ-Bild und Lage der Niveaulinien {g = d}
print("\n(4) Schlitz ζ = L(1+u) im τ-Bild (arg τ > 90° hiesse: neigt zum D-Bogen, < 90°: zum A-Bogen):")
for u in (1e-8, 1e-4, 1e-2, 0.05, 0.15):
    t = tau(L * (1 + mpf(u)))
    print(f"   u={u}: arg τ = {mp.nstr(arg(t) * 180 / pi, 6)}°, |τ| = {mp.nstr(abs(t), 4)}, g = {mp.nstr(g(L * (1 + mpf(u))), 4)}")
for d in (0.05, 0.1):
    dmin, dmax = mpf(10), mpf(0)
    for k in range(72):
        e = exp(2j * pi * k / 72)
        c0 = a0 / 2

        def tk(m, r):
            b = ((c0 - m) * conj(e)).real
            return -b + mp.sqrt(b ** 2 - (abs(c0 - m) ** 2 - r ** 2))

        r0 = max(tk(0, RHO), tk(a0, R0))
        try:
            x = findroot(lambda x: g(c0 + (r0 + x) * e) - d, mpf("0.1"))
        except Exception:
            continue
        dmin, dmax = min(dmin, x), max(dmax, x)
    print(f"   Niveaulinie g = {d}: radialer Abstand von ∂W (Strahlen von a0/2) zwischen {mp.nstr(dmin, 3)} und {mp.nstr(dmax, 3)}")
