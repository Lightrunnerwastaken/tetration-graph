"""Strang V — rigorose Geometrie für Satz K1 mit Linse (arb; f_K aus Vs (K2)-Datei über mk_kneser_arb), Auftrag Leitung 23:20, Punkt 1.

Normierte Variable z = x/ρ, F(z) = f_K(ρz)/ρ, e^{iθ} = L/ρ, μ = √L = F'(e^{iθ}). Linse z(τ) = e^{s(τ)}e^{i(θ+τ)}, s(τ) = S(1 − e^{−(τ/τ₁)^p}).
Geprüft (oberer Halbkreis; unterer per Spiegelung, f_K reell):
 (a) f(ℬ) ⊂ D: |F| < 1 auf Γin und auf der Linsenkurve (Maximumprinzip auf ℬ), Intervalle für τ ∈ [τ_min, π − θ];
     für τ ∈ (0, τ_min] lokal: F(z) = e^{iθ} + μΔ + R, |R| ≤ ½C₂|Δ|², C₂ = sup|F''| auf dem Segment [e^{iθ}, z] über die Darstellung
     F(z) = log(ρF(w))/ρ, w = e^{ρz}/ρ: F'' = ρw[F''(w)w + F'(w) − F'(w)²w/F(w)]/F(w), |w| ≤ 1 ⇒ C₂ ≤ ρ(S2 + S1 + S1²/m)/m, m = min|F(w)|.
 (b) f(Γout) trifft den Strahl S_L = {L(1+u), u > 0} nicht: arg F(e^{iφ}) < θ für φ ∈ [0, θ − τ_min] (Intervalle) und lokal für φ ∈ (θ − τ_min, θ);
     mit (a) und f(Γin) ⊂ D liegt S_L in der unbeschränkten Komponente von ℂ∖f(∂D) ⇒ S_L ∩ f(D̄) = ∅; Linse ∩ S_L = ∅ trivial (s(τ) > 0 nur für τ > 0).
 (c) M_γ = max_{D̄∪ℬ̄}|m| ≤ max über Γout ∪ Linsenkurve von |e^{ρz}/F'(z)| (Maximumprinzip; m = f'∘f = eˣ/f' analytisch auf D̄ ∪ ℬ̄).
Aufruf: python v32_geometrie_zert.py S [tau1] [p] [tau_min_grad]
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
TAU1 = arb(sys.argv[2] if len(sys.argv) > 2 else "12") * arb.pi() / 180
PEXP = arb(sys.argv[3] if len(sys.argv) > 3 else "1.3")
TMIN = arb(sys.argv[4] if len(sys.argv) > 4 else "0.2") * arb.pi() / 180
KN = KneserArb(str(HIER / "kneser_halbexp_b_zert_V.json"))
ctx.prec = 120
pi = arb.pi()
th, rho, mu = KN.th, KN.rho, KN.mu
eith = acb(th.cos(), th.sin())
S1, S2 = KN.S1, KN.S2
MCAP = arb(sys.argv[5] if len(sys.argv) > 5 else "2.4")      # Ziel-Schranke für |m| (erzwingt feinere Intervalle)


def s_of(t):
    return S * (1 - (-(t / TAU1) ** PEXP).exp())


def lens_ball(t_lo, t_hi):
    """Kugel, die die Linsenkurve über τ ∈ [t_lo, t_hi] enthält (Mittelpunkt + Radius aus Ableitungsschranke)."""
    tm = (t_lo + t_hi) / 2
    h = (t_hi - t_lo) / 2
    zc = s_of(tm).exp() * acb((th + tm).cos(), (th + tm).sin())
    # |dz/dτ| ≤ e^{S}(1 + max|s'|), s' ≤ S·p/τ₁·(τ/τ₁)^{p−1} e^{−(τ/τ₁)^p} ≤ S·p/τ₁ (grob)
    dz = S.exp() * (1 + S * PEXP / TAU1 * arb(max(1, float(((t_hi / TAU1) ** (PEXP - 1)).upper()))))
    r = (dz * h).upper()
    zc0 = acb(zc.real.mid(), zc.imag.mid())
    return acb(arb(zc0.real.mid(), r + 1e-30), arb(zc0.imag.mid(), r + 1e-30))


def zerlege(lo, hi, test, tiefe=0, maxtiefe=14):
    """Bisektion: test(lo, hi) → True (bestanden), sonst teilen."""
    if test(lo, hi):
        return 1, 0
    if tiefe >= maxtiefe:
        return 0, 1
    mid = (lo + hi) / 2
    a1, b1 = zerlege(lo, mid, test, tiefe + 1, maxtiefe)
    a2, b2 = zerlege(mid, hi, test, tiefe + 1, maxtiefe)
    return a1 + a2, b1 + b2


stats = {}
# (a) Γin: |F(e^{iφ})| < 1, φ = θ + τ
def test_gin(lo, hi):
    phi = arb(((th + lo + th + hi) / 2).mid(), ((hi - lo) / 2).upper() * (1 + 1e-9))
    try:
        v = KN.F_kreis(phi)
    except AssertionError:
        return False
    return bool(abs(v) < 1)


ok, fail = zerlege(TMIN, pi - th, test_gin)
stats["Gin"] = (ok, fail)
print(f"(a1) f(Γin) ⊂ D für τ ∈ [{float(TMIN * 180 / pi):.2f}°, {float((pi - th) * 180 / pi):.1f}°]: {ok} Intervalle bestanden, {fail} offen", flush=True)
M_lens = arb(0)


# (a) Linse: |F(z)| < 1 und gleich M_γ-Beitrag
def test_lens(lo, hi):
    global M_lens
    z = lens_ball(lo, hi)
    try:
        Fz, w, Fw = KN.F_gin(z)
        Fpz = KN.Fp_gin(z)
    except AssertionError:
        return False
    mval = abs((rho * z).exp() / Fpz)
    if not mval.is_finite():
        return False
    if not bool(abs(Fz) < 1) or not bool(mval < MCAP):
        return False
    M_lens = arb(max(M_lens.upper(), mval.upper()))
    return True


ok, fail = zerlege(TMIN, pi - th, test_lens)
stats["Linse"] = (ok, fail)
print(f"(a2) f(Linse) ⊂ D für τ ∈ [{float(TMIN * 180 / pi):.2f}°, …]: {ok} bestanden, {fail} offen; max|m| auf der Linse ≤ {float(M_lens):.4f}", flush=True)
# (a3) lokal τ ≤ τ_min: Δ = z − e^{iθ}; |e^{iθ} + μΔ| + ½C₂|Δ|² + ε < 1 ?
eps = arb(0)                             # F(e^{iθ}) = e^{iθ}, F'(e^{iθ}) = μ exakt (Satz p32a / (K3))
# C₂: w-Segment in |w − e^{iθ}| ≤ δ_w, m ≥ 1 − S1·δ_w (|F(e^{iθ})| = 1)
dz_max = TMIN * (1 + S)                 # |Δ| ≤ τ + s ≤ τ(1 + S) für τ ≤ τ_min (s ≤ S·τ/τ₁·… ≤ Sτ grob, da τ₁ > 1 rad? — ersetzt unten durch direkte Schranke)
dz_max = TMIN + s_of(TMIN) + arb("1e-30")
dw = rho * dz_max * arb(1).exp() ** (rho * dz_max)          # |w − e^{iθ}| ≤ |e^{ρΔ} − 1| ≤ ρ|Δ|e^{ρ|Δ|}
m_low = 1 - S1 * dw
C2 = rho * (S2 + S1 + S1 ** 2 / m_low) / m_low * (1 + dw)
lok_ok = True
nk = 400
schlechteste = arb(-10)
for k in range(nk):
    t_lo = TMIN * k / nk
    t_hi = TMIN * (k + 1) / nk
    # obere Schranke von |e^{iθ} + μΔ| über das Teilintervall: Δ(τ) = z(τ) − e^{iθ}
    tb = arb(((t_lo + t_hi) / 2).mid(), ((t_hi - t_lo) / 2).upper() * (1 + 1e-9))
    if bool(tb.lower() <= 0):
        tb = arb(((t_lo + t_hi) / 2).mid(), ((t_hi - t_lo) / 2).upper())
    z = s_of(tb).exp() * acb((th + tb).cos(), (th + tb).sin())
    D = z - eith
    lin = abs(eith + mu * D)
    val = lin + C2 * abs(D) ** 2 / 2 + eps * (1 + abs(D))
    # Für τ → 0 ist val → 1; prüfe val < 1 auf τ ∈ [t_lo, t_hi] mit t_lo > 0 strikt; das Intervall [0, τ_min/nk] gesondert (Steigung)
    if k == 0:
        # Ableitung nach τ bei 0: d|e^{iθ}+μΔ|/dτ = Re(e^{−iθ}μ·i) = −Im μ < 0; Restglied O(τ²), lineare Abnahme dominiert für τ ≤ τ_min/nk
        steig = -mu.imag + C2 * (1 + S) * (t_hi) + S * PEXP / TAU1 * (t_hi / TAU1) ** (PEXP - 1) * abs(mu)
        if not bool(steig < 0):
            lok_ok = False
        continue
    schlechteste = arb(max(schlechteste.upper(), (val - 1).upper()))
    if not bool(val < 1):
        lok_ok = False
print(f"(a3) lokal τ ∈ (0, τ_min]: C₂ ≤ {float(C2.upper()):.3f}, m ≥ {float(m_low.lower()):.4f}; max(|F| − 1) ≤ {float(schlechteste):.3e} — "
      f"{'bestanden' if lok_ok else 'NICHT bestanden'}", flush=True)
# (b) arg F(e^{iφ}) < θ für φ ∈ [0, θ − τ_min]; dazu M_γ auf Γout
M_out = arb(0)


def test_out(lo, hi):
    global M_out
    phi = arb(((lo + hi) / 2).mid(), ((hi - lo) / 2).upper() * (1 + 1e-9))
    try:
        v = KN.F_kreis(phi)
        vp = KN.Fp_kreis(phi)
    except AssertionError:
        return False
    z = acb(phi.cos(), phi.sin())
    mval = abs((rho * z).exp() / vp)
    if not mval.is_finite() or not bool(mval < MCAP):
        return False
    # arg v < θ ⇔ Im(v·e^{−iθ}) < 0 (v in der Nähe von e^{iθ}, Re(v e^{−iθ}) > 0) — allgemein: arg v ∈ (−θ, θ)
    q = v * acb(th.cos(), -th.sin())          # Im q < 0 ⇒ v liegt nicht auf dem Strahl arg = θ
    if not bool(q.imag < 0):
        return False
    M_out = arb(max(M_out.upper(), mval.upper()))
    return True


ok, fail = zerlege(arb(0), th - TMIN, test_out)
print(f"(b) arg f(e^{{iφ}}) < θ für φ ∈ [0, θ − τ_min]: {ok} bestanden, {fail} offen; max|m| auf Γout ≤ {float(M_out):.4f}", flush=True)
# (b) lokal: F(e^{i(θ−τ)}) − e^{iθ} = μΔ + R, Δ = e^{i(θ−τ)} − e^{iθ}; Im((μΔ)e^{−iθ}) = Im(μ(e^{−iτ} − 1)) ≈ −τ Re μ < 0
lok_b = True
# erstes Teilintervall (0, τ_min/nk]: Im q ≤ −τ·Re μ·(sin τ/τ) + |μ|τ²/2 + C₂τ²/2 < 0, da Re μ > (|μ| + C₂)·τ_min/(2nk)
if not bool(mu.real > (abs(mu) + C2) * TMIN / (2 * nk) + 0):
    lok_b = False
for k in range(2, nk + 1):
    tb = arb((TMIN * (k - arb("0.5")) / nk).mid(), (TMIN / (2 * nk)).upper() * (1 + 1e-9))
    D = acb((th - tb).cos(), (th - tb).sin()) - eith
    q = (mu * D) * acb(th.cos(), -th.sin())
    if not bool(q.imag + C2 * abs(D) ** 2 / 2 < 0):
        lok_b = False
print(f"(b) lokal φ ∈ (θ − τ_min, θ): {'bestanden' if lok_b else 'NICHT bestanden'}", flush=True)
print(f"(c) M_γ ≤ max({float(M_out):.4f}, {float(M_lens):.4f}) = {max(float(M_out), float(M_lens)):.4f}", flush=True)
