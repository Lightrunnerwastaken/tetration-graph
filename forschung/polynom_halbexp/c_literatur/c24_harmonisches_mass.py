"""Strang C, Runde 3, c24: harmonisches Mass von E = Ĉ∖W̄ — Poisson–Jensen-Variante von Satz 12.2.

log R_N(ζ) ≤ ∫_{∂W} log R_N dω_ζ − (N+2) g_E(ζ).  Ist log R_N(η) ≤ κ log N + log C + N Λ(η) auf ∂W, so ist
log R_N(ζ) ≤ κ log N + log C + N[ω_ζ(Λ) − g_E(ζ)] − 2 g_E(ζ),  ω_ζ(Λ) := ∫ Λ dω_ζ.
Hier: (1) Randparameter τ ∈ ℝ (Ecke L ↦ 0, Ecke L̄ ↦ ∞) der Punkte von c23; (2) für Λ = λ·1_{D-Bogen jenseits t ≥ t1} das Verhältnis
ω_ζ(Λ)/g_E(ζ) längs des Schlitzes ζ = L(1+u) und auf der Kontur ∂(W+0.1); (3) kritisches λ, unter dem die Konsistenzschranke bleibt.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, log, arg, conj, atan

mp.dps = 40
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


def gE(z):
    t = tau(z)
    return log(abs((t - conj(TAU0)) / (t - TAU0)))


def omega_intervall(z, a, b):
    """harmonisches Mass (in E, Pol ∞ beachtet über die Halbebene) von τ ∈ (a,b) ⊂ ℝ, gesehen von z."""
    t = tau(z)
    x, y = t.real, t.imag
    fa = atan((a - x) / y) if a != -mp.inf else -pi / 2
    fb = atan((b - x) / y) if b != mp.inf else pi / 2
    return (fb - fa) / pi


aL = arg(L)
print(f"α = {mp.nstr(ALPHA, 7)}, τ0 = τ(∞) = {mp.nstr(TAU0, 6)}")
# (1) Randparameter der D-Bogen-Punkte
tauD = {}
for t in [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.0]:
    z = RHO * exp(1j * (aL + t * (pi - aL))) * (1 + mpf(10) ** -25)
    tauD[t] = tau(z).real
print("τ der D-Bogen-Punkte (t: τ):", " ".join(f"{t}:{mp.nstr(v, 4)}" for t, v in tauD.items()))
# Seite des D-Bogens in τ: Vorzeichen
seite = 1 if tauD[1.0] > 0 else -1
# (2) Verhältnis ω(D-Bogen jenseits t1)/g_E längs des Schlitzes
for t1 in (0.02, 0.1, 0.2):
    s1 = abs(tauD[t1])
    zeile = []
    for u in (1e-4, 1e-3, 1e-2, 0.05, 0.2):
        z = L * (1 + u)
        om = omega_intervall(z, s1, mp.inf) if seite > 0 else omega_intervall(z, -mp.inf, -s1)
        zeile.append(f"u={u}: {mp.nstr(om / gE(z), 4)}")
    print(f"ω_ζ(D-Bogen, t ≥ {t1}) / g_E(ζ) längs Schlitz: " + "  ".join(zeile))
# ganzer D-Bogen
zeile = []
for u in (1e-4, 1e-3, 1e-2, 0.05, 0.2):
    z = L * (1 + u)
    om = omega_intervall(z, 0, mp.inf) if seite > 0 else omega_intervall(z, -mp.inf, 0)
    zeile.append(f"u={u}: ω={mp.nstr(om, 4)}, ω/g_E={mp.nstr(om / gE(z), 4)}")
print("ganzer D-Bogen: " + "  ".join(zeile))
# (3) Kontur ∂(W+d): min g_E und max ω(D-Bogen)
for d in (0.05, 0.1, 0.2):
    mn, mx = mpf(10), mpf(0)
    for k in range(144):
        phi = 2 * pi * k / 144
        e = exp(1j * phi)
        c0 = a0 / 2
        def tk(m, r):
            b = ((c0 - m) * e.conjugate()).real
            return -b + mp.sqrt(b ** 2 - (abs(c0 - m) ** 2 - r ** 2))
        z = c0 + (max(tk(0, RHO), tk(a0, R0)) + d) * e
        if abs(z - L) < 0.25 or abs(z - Lq) < 0.25:
            continue
        g = gE(z)
        om = omega_intervall(z, 0, mp.inf) if seite > 0 else omega_intervall(z, -mp.inf, 0)
        mn = min(mn, g / max(om, mpf(10) ** -30))
    print(f"Kontur ∂(W+{d}) ohne |ζ−L|<0.25: min g_E/ω(D-Bogen) = {mp.nstr(mn, 4)}  (⇒ λ < diesem Wert genügt dort)")
