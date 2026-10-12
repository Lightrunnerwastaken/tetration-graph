"""Strang V — Green-Konstanten von Satz V2 streng (Auftrag Leitung 02:05, Punkt 2), Teil 1: Vergleichsgebiet „Lune“.

Idee (Vergleichsprinzip): Für K_γ ⊂ K* gilt g_{K_γ} ≥ g_{K*} ausserhalb K*. Als K* dient die Vereinigung zweier Kreisscheiben D₁ ∪ D₂, deren Ränder
beide durch L und L̄ gehen (Mittelpunkte c₁ < 0 < c₂ reell). Ihr Komplement ist unter T(z) = (z − L)/(z − L̄) ein Sektor der Öffnung απ, also
g_{K*} explizit: ζ = (e^{−i·lo}T(z))^{1/α}, g_{K*}(z) = log|(ζ − ζ̄_∞)/(ζ − ζ_∞)|, ζ_∞ = Bild von ∞ (T(∞) = 1).
Geprüft (arb, f_K aus Vs (K2)-Datei, Linse Ks Form):
 (i)  K_γ ⊂ D₁ ∪ D₂: Linsenkurve ⊂ D₁ (analytisch, Intervalle); Γin ⊂ D₁ und Γout ⊂ D₂ (Kreisbüschel durch L, L̄: c₁ < 0 < c₂);
      f(Γout) ⊂ D₂ (Intervalle über F_kreis, nahe L lokal: F = e^{iθ} + μΔ + R); f(Γin), f(ℬ) ⊂ D ⊂ D₁ ∪ D₂ (v32). D₁ ∪ D₂ einfach zusammenhängend
      ⇒ das von f(∂D) umschlossene f(D̄) liegt darin.
 (ii) Schranken längs des Schlitzes: g_γ(L(1+u)) ≥ g_{K*}(L(1+u)) — als Tabelle (arb) für u ∈ [1e-4, 0.4].
 (iii) Niveaulinie Γ*_d = {g_{K*} = d}: dort g_γ ≥ d (streng); Länge, Abstand zu D̄, Schnitt u₀ mit dem Schlitz (arb aus der expliziten Parametrisierung).
Aufruf: python v39_lune_zert.py S c1 c2 [d]
"""
import json
import math
import sys
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
from mk_kneser_arb import KneserArb  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
S = arb(sys.argv[1])
c1, c2 = arb(sys.argv[2]), arb(sys.argv[3])
DS = [arb(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "0.08,0.09,0.1").split(",")]
import os
TAU1 = arb(os.environ.get("V_TAU1", "12")) * arb.pi() / 180
PEXP = arb(os.environ.get("V_P", "1.3"))
KN = KneserArb(str(HIER / "kneser_halbexp_b_zert_V.json"))
ctx.prec = 120
pi = arb.pi()
th, rho, mu = KN.th, KN.rho, KN.mu
L = rho * acb(th.cos(), th.sin())
Lb = L.conjugate()
r1, r2 = abs(c1 - L), abs(c2 - L)


def s_of(t):
    return S * (1 - (-(t / TAU1) ** PEXP).exp())


def T(z):
    return (z - L) / (z - Lb)


# Sektorgrenzen: Bild eines Kreispunkts (Gerade durch 0)
def gwinkel(c, r):
    p = c + r * acb(arb("0.3").cos(), arb("0.3").sin())
    a = T(p).arg()
    while bool(a < 0):
        a = a + pi
    while bool(a >= pi):
        a = a - pi
    return a


g1, g2 = gwinkel(c1, r1), gwinkel(c2, r2)
lo, hi = g2 - pi, g1
alpha = (hi - lo) / pi
print(f"[S = {float(S.mid())}] D₁: c = {float(c1.mid())}, r = {float(r1.mid()):.5f}; D₂: c = {float(c2.mid())}, r = {float(r2.mid()):.5f}; "
      f"Sektor ({float((lo * 180 / pi).mid()):.3f}°, {float((hi * 180 / pi).mid()):.3f}°), α = {float(alpha.mid()):.5f}", flush=True)
rot = acb((-lo).cos(), (-lo).sin())


def zeta(z):
    w = T(z) * rot
    return (w.log() / alpha).exp()


zinf = ((rot).log() / alpha).exp()


def g_lune(z):
    zz = zeta(z)
    return ((zz - zinf.conjugate()) / (zz - zinf)).abs_lower().log() if False else abs((zz - zinf.conjugate()) / (zz - zinf)).log()


def zerlege(lo_, hi_, test, tiefe=0, maxt=16):
    if test(lo_, hi_):
        return 1, 0
    if tiefe >= maxt:
        return 0, 1
    m = (lo_ + hi_) / 2
    a1, b1 = zerlege(lo_, m, test, tiefe + 1, maxt)
    a2, b2 = zerlege(m, hi_, test, tiefe + 1, maxt)
    return a1 + a2, b1 + b2


# (i) Linsenkurve ⊂ D₁: |ρe^{s}e^{iφ} − c₁| < r₁, φ = θ + τ
def test_lin(t_lo, t_hi):
    tb = arb(((t_lo + t_hi) / 2).mid(), ((t_hi - t_lo) / 2).upper() * (1 + 1e-9))
    z = rho * s_of(tb).exp() * acb((th + tb).cos(), (th + tb).sin())
    return bool(abs(z - c1) < r1)


TMIN = arb("1e-3")
ok, fail = zerlege(TMIN, pi - th, test_lin)
# nahe L: |z − c₁|² − r₁² = 2Re((L − c₁)·conj(Δ)) + |Δ|², Δ = z − L ≈ ρe^{iθ}(iτ + s); Vorzeichen über die Ableitung
D_dir = acb(0, 1) * L                    # dz/dτ bei τ = 0 (s'(0) = 0 für p > 1)
if bool(PEXP == 1):
    # Keil (p = 1): s'(0) = S/τ₁ ≠ 0; streng auf [0, TMIN]: h(τ) = |z − c₁|² − r₁² ≤ h'(0)τ + M₂τ²/2 mit M₂ ≥ sup|h''| (Intervalle)
    D_dir = L * acb(S / TAU1, 1)
    h1 = 2 * ((L - c1).conjugate() * D_dir).real
    tb = arb(TMIN / 2, TMIN / 2)
    zb = rho * s_of(tb).exp() * acb((th + tb).cos(), (th + tb).sin())
    sp = S / TAU1 * (-tb / TAU1).exp()
    spp = -S / TAU1 ** 2 * (-tb / TAU1).exp()
    z1 = zb * acb(sp, 1)
    z2 = zb * (acb(sp, 1) ** 2 + spp)
    M2h = (2 * abs(z1) ** 2 + 2 * abs(zb - c1) * abs(z2)).upper()
    lok = bool(h1 < 0) and bool(h1 + M2h * TMIN / 2 < 0)
    print(f"    Keil: h'(0) = {h1.str(5)}, sup|h''| ≤ {arb(M2h).str(4)} auf [0, {TMIN.str(2)}] ⇒ h < 0 dort: {lok}", flush=True)
else:
    lok = bool((2 * ((L - c1) * D_dir.conjugate()).real) < 0)
print(f"(i-a) Linse ⊂ D₁: {ok} Intervalle für τ ∈ [1e-3, π − θ], {fail} offen; nahe L Ableitung von |z − c₁|² < 0: {lok}", flush=True)


# (i) f(Γout) ⊂ D₂: |ρF(e^{iφ}) − c₂| < r₂ für φ ∈ [0, θ − τ_min]
def test_fout(p_lo, p_hi):
    phi = arb(((p_lo + p_hi) / 2).mid(), ((p_hi - p_lo) / 2).upper() * (1 + 1e-9))
    try:
        v = KN.F_kreis(phi) * rho
    except AssertionError:
        return False
    return bool(abs(v - c2) < r2)


ok2, fail2 = zerlege(arb(0), th - TMIN, test_fout)
# nahe L: f(e^{i(θ−τ)}) − L = μρ(e^{i(θ−τ)} − e^{iθ}) + R; Richtung −iμL; |f − c₂|² − r₂² ≈ 2Re((L − c₂)conj(−iμLτ)) < 0 ?
lok2 = bool((2 * ((L - c2) * (acb(0, -1) * mu * L).conjugate()).real) < 0)
print(f"(i-b) f(Γout) ⊂ D₂: {ok2} Intervalle für φ ∈ [0, θ − 1e-3], {fail2} offen; nahe L lineare Richtung nach innen: {lok2}", flush=True)
# (ii) Schlitz
tab = []
for k in range(241):
    u = arb(10) ** (-4 + arb(k) * arb("3.6") / 240)
    gl = g_lune(L * (1 + u))
    tab.append([float(u.mid()), float(gl.lower())])
print("(ii) g_γ(L(1+u)) ≥ g_{K*}: " + ", ".join(f"u = {u:.1e}: {g:.4e}" for u, g in tab[::48]) + f"; monoton steigend: {all(tab[i][1] < tab[i + 1][1] for i in range(len(tab) - 1))}", flush=True)
# (iii) Niveaulinien Γ*_d: in der ζ-Ebene |ζ − ζ̄∞| = e^d|ζ − ζ∞| ist ein Kreis (Apollonius); Parametrisierung, Rücktransformation
niv = {}
for d in DS:
    q = d.exp()
    # Apollonius-Kreis: Punkte ζ mit |ζ − a| = q|ζ − b|, a = ζ̄∞, b = ζ∞ ⇒ Mittelpunkt m = (a − q²b)/(1 − q²), Radius R = q|a − b|/|1 − q²|
    a_, b_ = zinf.conjugate(), zinf
    m = (a_ - q * q * b_) / (1 - q * q)
    R = q * abs(a_ - b_) / abs(1 - q * q)
    npkt = 6000
    pts = []
    laenge = arb(0)
    prev = None
    u0 = None
    minabs = arb(100)
    maxabs = arb(0)
    for j in range(npkt + 1):
        t = 2 * pi * j / npkt
        zz = m + R * acb(t.cos(), t.sin())
        if not bool(zz.imag > 0):           # nur der Teil in der oberen Halbebene (Sektorbild)
            prev = None
            continue
        w = (zz.log() * alpha).exp() / rot
        z = (L - Lb * w) / (1 - w)            # T⁻¹
        pts.append([float(z.real.mid()), float(z.imag.mid())])
        if prev is not None:
            laenge += abs(z - prev)
        prev = z
        minabs = arb(min(minabs.upper(), abs(z).lower()))
        maxabs = arb(max(maxabs.upper(), abs(z).upper()))
    # Schnitt mit dem Schlitz: u mit g_{K*}(L(1+u)) = d (Bisektion)
    ulo, uhi = arb("1e-4"), arb("2")
    for _ in range(60):
        um = (ulo + uhi) / 2
        if bool(g_lune(L * (1 + um)) < d):
            ulo = um
        else:
            uhi = um
    niv[str(float(d.mid()))] = dict(laenge=float(laenge.upper()) * 1.001, min_abs=float(minabs.lower()), max_abs=float(maxabs.upper()),
                                    u0=float(uhi.upper()), punkte=pts)
    print(f"(iii) Γ*_d, d = {float(d.mid())}: Länge ≤ {float(laenge.upper()) * 1.001:.3f}, min|x| = {float(minabs.lower()):.4f}, "
          f"max|x| = {float(maxabs.upper()):.3f}, Schlitz-Schnitt u₀ = {float(uhi.upper()):.4f}", flush=True)
json.dump(dict(S=float(S.mid()), c1=float(c1.mid()), c2=float(c2.mid()), alpha=float(alpha.mid()), lo=float(lo.mid()), schlitz=tab, niveau=niv,
               linse_ok=fail == 0 and lok, fout_ok=fail2 == 0 and lok2), open(HIER / f"v39_lune_S{float(S.mid())}_t{os.environ.get('V_TAU1', '12')}_p{os.environ.get('V_P', '1.3')}.json", "w"))
