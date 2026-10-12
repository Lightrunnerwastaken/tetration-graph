"""Strang V — Ks E2 (BERICHT §24.5): Eckbedingung (EB_N) an den echten p_N messen.

(EB_N): |p_N − f_K| ≤ η_N auf 𝒲_N = {|x − L| ≤ A₀r_N} ∪ ⋃_{r_N ≤ r ≤ ε}(Ruhefenster bei Radius r), η_N ≪ r_N, η_N·ln N ≪ 1; r_N = ρN^{−Re a}.
f_K nahe L (auch ausserhalb D̄, Schlitzebene mit Schnitt S_L) aus Ps oberer Struktur wie v44: W = P⁻¹(x), arg W ∈ (θ + 4π, θ + 6π) (Zweig k₀ = 2,
in v44 über Ps slog-Krawczyk bestimmt), z + θ_K(z) = (ln|W| + i·arg W)/L, f_K(x) = P(√L·W·e^{L(θ_K(z+½) − θ_K(z))}). p_N: Strang A (N ≤ 1000, 60 Stellen),
V v33 (N = 1500, 2000). Ruhefenster: physikalischer Winkel θ_phys = u_A − arg Λ₀ − Im(a)·ln N (Grad), Λ₀ = −λ_Aρ/L², λ_A = −0.203297 − 1.406007i,
u_A ∈ (19°, 100°) ∪ (−120°, −61°) (Ks §24.5).
Aufruf: python v49_E2_EB.py [N-Liste]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "pv_kopie"))
from kc_kneser import Kneser  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "200,400,600,800,1000,1500,2000").split(",")]
PREC = 260
KN = Kneser(PREC, laut=False)
R = KN.R
ctx.prec = PREC
L = R.L
c = KN.c
pi = arb.pi()
TPI = 2 * pi * acb(0, 1)
th = L.imag
sqL = (L / 2).exp()
codd = [c[n] if n % 2 == 1 else acb(0) for n in range(len(c))]
Lc = complex(float(L.real.mid()), float(L.imag.mid()))
rho = abs(Lc)
a_exp = 1 + 1j * Lc / (4 * math.pi)
lamA = complex(-0.203297, -1.406007)
Lam0 = -lamA * rho / Lc ** 2
argLam0 = math.degrees(cmath.phase(Lam0))


def mid(z):
    return acb(z.real.mid(), z.imag.mid())


def fK(x):
    """f_K(x) nahe L (Schlitzebene), x acb."""
    w = mid(x - L)
    for _ in range(60):
        f, df = R.P_dP_ball(w)
        dw = mid((f - x) / df)
        w = mid(w - dw)
        if dw.abs_upper() < arb(2) ** (-PREC + 20):
            break
    A = w.arg()
    while bool(A <= th + 4 * pi):
        A += 2 * pi
    while bool(A > th + 6 * pi):
        A -= 2 * pi
    zeta = acb(abs(w).log(), A) / L
    z = zeta - c[0]
    for _ in range(100):
        zn = mid(zeta - R._horner(c, (TPI * z).exp()))
        if (zn - z).abs_upper() < arb(2) ** (-PREC + 20):
            z = zn
            break
        z = zn
    assert float(z.imag.mid()) > 0.36, "Im z < y₁"
    d = -2 * R._horner(codd, (TPI * z).exp())
    return R.P_ball(sqL * w * (L * d).exp())


def lade_pN(N):
    ctx.prec = 900                      # Horner nahe |x| = ρ braucht bei N ≥ 1500 > 260 Bit (gemessen: 260 Bit ⇒ Radius 2e-6, 600 Bit ⇒ 5e-21)
    try:
        return _lade(N)
    finally:
        ctx.prec = PREC


def _lade(N):
    if N <= 1000:
        d = json.loads((HIER.parent.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"][str(N)]
    else:
        d = json.loads((HIER / "sg" / f"v33_pN_{N}.json").read_text(encoding="utf-8"))["p"]
    return [acb(arb(s)) for s in d]


def horner(cf, x):
    alt_ = ctx.prec
    ctx.prec = 900
    s = acb(0)
    x = acb(x.real.mid(), x.imag.mid())     # Mittelpunkt: acb-Radien wachsen je Horner-Schritt um (|Re x| + |Im x|)/|x| ≈ 1.2 (nahe L), bei N = 2000 ⇒ 1e158
    for a in reversed(cf):
        s = s * x + a
    ctx.prec = alt_
    return s


# Kontrolle von f_K gegen die zertifizierte Taylor-Reihe in D (P-V-Koeffizienten)
zk = json.loads((HIER.parent.parent / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))["a"]
ak = [acb(arb(x[0] if isinstance(x, list) else x)) for x in zk]
for xx in (L * arb("0.93"), L + 0.1 * acb(math.cos(th + 2.5), math.sin(th + 2.5))):
    print(f"Kontrolle f_K: obere Struktur {mid(fK(xx)).str(15)} | Taylor-Reihe (641 Terme) {mid(horner(ak, xx)).str(15)}", flush=True)
t0 = time.time()
erg = {}
for N in NS:
    try:
        pN = lade_pN(N)
    except (FileNotFoundError, KeyError):
        print(f"N = {N}: keine Daten", flush=True)
        continue
    rN = rho * N ** (-a_exp.real)
    # (1) Kegelscheibe |x − L| ≤ 16 r_N, Polargitter; Schlitzsektor ±3° um θ getrennt
    sup_aus, sup_schlitz, sup_innen = 0.0, 0.0, 0.0
    thg = math.degrees(float(th.mid()))
    sek = {}                                       # (Sektor, Radiusband) → max
    for i in range(14):
        r = 16 * rN * 10 ** (-3 + 3 * i / 13)
        band = "≤4r_N" if r <= 4 * rN * 1.0001 else "4–16r_N"
        for k in range(120):
            ang = 3.0 * k + 1.5
            x = L + arb(r) * acb(math.cos(math.radians(ang)), math.sin(math.radians(ang)))
            e = float(abs(horner(pN, x) - fK(x)).mid())
            dang = (ang - thg) % 360                # 0 = Schlitzrichtung, 90 = Γin-Tangente, 270 = Γout-Tangente
            if dang < 3.0 or dang > 357.0:
                sup_schlitz = max(sup_schlitz, e)
            else:
                sup_aus = max(sup_aus, e)
            sektor = "Linsenseite" if dang < 90 else ("D" if dang <= 270 else "Γout-Seite")
            sek[(sektor, band)] = max(sek.get((sektor, band), 0.0), e)
            if abs(complex(float(x.real.mid()), float(x.imag.mid()))) <= rho:
                sup_innen = max(sup_innen, e)
    print("   Sektoren (aussen Linsenseite = zwischen S_L und Γin, aussen Γout-Seite = zwischen Γout und S_L): " +
          ", ".join(f"{a} {b}: {v:.2e}" for (a, b), v in sorted(sek.items())), flush=True)
    # (2) Ruhefenster-Strahlen
    strahlen = {}
    for uA in (25, 45, 65, 85, -115, -95, -75):
        phi = uA - argLam0 - math.degrees(a_exp.imag) * math.log(N)
        mx, mx16, prof = 0.0, 0.0, []
        for i in range(25):
            r = rN * (0.05 / rN) ** (i / 24)
            x = L + arb(r) * acb(math.cos(math.radians(phi)), math.sin(math.radians(phi)))
            e = float(abs(horner(pN, x) - fK(x)).mid())
            mx = max(mx, e)
            if r <= 16 * rN:
                mx16 = max(mx16, e)
            prof.append((r, e))
        strahlen[uA] = (phi % 360, mx)
        print(f"   Strahl u_A = {uA}° (φ = {phi % 360:.1f}°): max auf [r_N, 16r_N] {mx16:.2e}; Profil " +
              " ".join(f"r={r:.1e}:{e:.1e}" for r, e in prof[::4]), flush=True)
    erg[N] = dict(rN=rN, kegel_ohne_schlitz=sup_aus, schlitzsektor=sup_schlitz, kegel_in_D=sup_innen, strahlen=strahlen)
    print(f"N = {N} ({time.time() - t0:.0f}s): r_N = {rN:.3e}, N^-4.8 = {N ** -4.8:.2e} | Kegelscheibe |x−L| ≤ 16r_N: sup|p_N − f_K| = {sup_aus:.2e} "
          f"(ohne ±3° um S_L), {sup_schlitz:.2e} (Schlitzsektor), {sup_innen:.2e} (Teil in D̄) | Ruhefenster r ∈ [r_N, 0.05]: " +
          ", ".join(f"u_A = {u}° (φ = {p:.1f}°): {m:.2e}" for u, (p, m) in strahlen.items()), flush=True)
json.dump({str(k): v for k, v in erg.items()}, open(HIER / "v49_E2.json", "w", encoding="utf-8"), ensure_ascii=False, indent=0)
