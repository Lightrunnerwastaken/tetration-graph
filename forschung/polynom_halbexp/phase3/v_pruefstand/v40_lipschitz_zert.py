"""Strang V — Green-Konstanten von Satz V2 streng, Teil 2: Lipschitz-Schranke für g_γ nahe ∂K_γ über innere Kreisscheiben.

Für η mit dist(η, K_γ) = δ und nächstem Randpunkt ζ: liegt eine Scheibe B vom Radius r₀ in K_γ und berührt ∂K_γ in ζ (innere Normale), so ist
g_γ(η) ≤ g_B(η) = log(1 + δ/r₀) ≤ δ/r₀ (Vergleichsprinzip). Also Lip ≤ 1/r₀.
r₀ aus Krümmungsschranken (Blaschke, rollende Scheibe): ∂K_γ = Linsenkurve ∪ f(Γout), stückweise C², in L, L̄ einspringende Ecke (Innenwinkel 218°),
dort keine Einschränkung. r₀ = min(1/κ_max⁺, halbe Engstelle).
  Linse r(φ) = ρe^{s}: κ = (1 + s'² − s'')/(ρe^{s}(1 + s'²)^{3/2}) (Intervalle in τ, s', s'' explizit).
  f(Γout): κ = Re(1 + z f''(z)/f'(z))/(|z||f'(z)|), z = ρe^{iφ}; f', f'' mit Horner + Koeffizienten-/Schwanzradien aus Vs (K2)-Datei (arb).
Engstelle: kleinster Abstand zwischen Linsenkurve und f(Γout) ausserhalb einer Umgebung von L, L̄ (Gleitkomma-Gitter plus Lipschitz-Rand).
Aufruf: python v40_lipschitz_zert.py S
"""
import math
import sys
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
from mk_kneser_arb import KneserArb  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
S = arb(sys.argv[1])
import os
TAU1 = arb(os.environ.get("V_TAU1", "12")) * arb.pi() / 180
PEXP = arb(os.environ.get("V_P", "1.3"))
KN = KneserArb(str(HIER / "kneser_halbexp_b_zert_V.json"))
ctx.prec = 120
pi = arb.pi()
th, rho = KN.th, KN.rho
Ke = KN.K_eval
b = KN.b
# Koeffizienten für F'' (normiert): F''(z) = (1/ρ)Σ k(k−1)b_k z^{k−2}
bk2 = [acb(arb((arb(k) * (k - 1) * b[k]).mid())) for k in range(2, Ke + 1)]
Rc2 = arb(sum((arb(k) * (k - 1) * b[k].rad() for k in range(2, Ke + 1)), arb(0)).upper())
R2 = arb(0)
for k in range(Ke + 1, KN.Kd + 1):
    R2 += abs(b[k]) * arb(k) ** 2
R2 += KN.C * arb(KN.Kd) ** (3 - KN.p) / (KN.p - 3)
R2 = arb(R2.upper())


def F2_kreis(phi):
    z0, r = KN._kreispunkt(phi)
    s = KN._horner(bk2, z0)
    e = Rc2 + R2 + r * KN.S3
    return (s + acb(arb(0, e.upper()), arb(0, e.upper()))) / KN.rho


def s_ab(t):
    """s, s', s'' für s = S(1 − e^{−(t/τ₁)^p})."""
    if bool(PEXP == 1):                       # Keil: glatt bis t = 0
        e = (-t / TAU1).exp()
        return S * (1 - e), S / TAU1 * e, -S / TAU1 ** 2 * e
    x = (t / TAU1) ** PEXP
    e = (-x).exp()
    s = S * (1 - e)
    s1 = S * e * PEXP * x / t
    s2 = S * e * (PEXP * (PEXP - 1) * x / t ** 2 - (PEXP * x / t) ** 2)
    return s, s1, s2


def zerlege(lo_, hi_, f, tiefe=0, maxt=14):
    """Obere Schranke von f über [lo, hi] durch Bisektion bis Auflösung; gibt max(obere Schranken) zurück."""
    tb = arb(((lo_ + hi_) / 2).mid(), ((hi_ - lo_) / 2).upper() * (1 + 1e-9))
    v = f(tb)
    if tiefe >= maxt or (v.rad() < 0.01 * abs(v.mid()) + 1e-6):
        return v.upper()
    m = (lo_ + hi_) / 2
    return max(zerlege(lo_, m, f, tiefe + 1, maxt), zerlege(m, hi_, f, tiefe + 1, maxt))


def kappa_lin(t):
    s, s1, s2 = s_ab(t)
    return (1 + s1 ** 2 - s2) / (rho * s.exp() * (1 + s1 ** 2) ** arb("1.5"))


if bool(PEXP == 1):
    km_lin = float(zerlege(arb(0), pi - th, kappa_lin))      # Keil: Intervalle bis τ = 0
else:
    km_lin = zerlege(arb("0.01") * pi / 180, pi - th, kappa_lin)
    # τ < 0.01°: s'' > 0 groß (p < 2) ⇒ κ < 1/ρ dort (Zähler 1 + s'² − s'' < 1 + s'²), Schranke 1/ρ·(1 + s'²) ≤ 1.01/ρ
    km_lin = max(float(km_lin), float((1.01 / rho).upper()))


def kappa_out(phi):
    z = rho * acb(phi.cos(), phi.sin())
    f1 = KN.Fp_kreis(phi)                     # f'(z) = F'(z/ρ)
    f2 = F2_kreis(phi) / rho
    return (1 + z * f2 / f1).real / (abs(z) * abs(f1))


km_out = float(zerlege(arb(0), th, kappa_out))          # bis einschliesslich φ = θ (Reihe konvergiert auf |z| = 1 dank (K2))
kmax = max(km_lin, km_out)
print(f"[S = {float(S.mid())}] κ_max Linse ≤ {km_lin:.4f}, κ_max f(Γout) ≤ {km_out:.4f}  ⇒  r₀ ≥ {1 / kmax:.4f} (x-Einheiten), "
      f"Lip(g_γ) ≤ 1/r₀ = {kmax:.4f}", flush=True)
