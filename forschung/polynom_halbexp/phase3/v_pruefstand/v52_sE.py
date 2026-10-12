"""Strang V — Rate s_E von dist_E(f_K, P_N) an der Spiralecke (Ks Satz K19, k_kern §30.2; Auftrag Leitung 10:55).

E := (D̄ ∖ U_{r₀}(L, L̄)) ∪ int𝒞 ∪ f_K(int𝒞), 𝒞 geschlossen um 0 durch L, L̄: Γout-Bogen |t| = ρ (bis Abstand r₁ von L) — Verbindungsbogen |t − L| = r₁ — Σ₂ nach L —
Σ₁ aus L — Verbindungsbogen — Linse |t| = ρ(1 + s_L) über Γin bis −ρ(1 + s_L), unten gespiegelt. Spiralen θ(r) = u_A + 71.46° + 1.623°·ln(r/ρ) (Ks §28.2).
Nahe L ist E der Spiralsektor von Σ₁ über die D-Seite bis f(Σ₂) (Drehung um arg μ = 38.3°): Aussenwinkel π·α_E = (u_A1 − u_A2 − 38.3°) (fest, dreht mit),
Vorhersage s_E = 5.447·α_E/(1 + T²), T = 0.02833 (Begradigung w = (z − L)^{1/(1+iT)}).
f_K auf ∂E: in D̄ Taylor (A, 4000 Koeff.), nahe L Ps obere Struktur (Zweig k₀ = 2, wie v44/v49), auf der Linse fern von L f(t) = Log f(eᵗ) (stetiger Zweig),
auf f(𝒞) exakt f_K(f(t)) = eᵗ. Bestapproximation diskret: Vandermonde mit Arnoldi auf den Randpunkten, Lawson-Iteration (Minimax), N = 20 … 400.
Aufruf: python v52_sE.py [uA1] [uA2] [r0]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "pv_kopie"))
sys.stdout.reconfigure(encoding="utf-8")
UA1 = float(sys.argv[1]) if len(sys.argv) > 1 else 80.0
UA2 = float(sys.argv[2]) if len(sys.argv) > 2 else -90.0
R0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.02
R1 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
SL, RMIN = 0.03, 1e-4
L = complex(0.31813150520476413, 1.3372357014306895)
rho, th = abs(L), cmath.phase(L)
TS = math.radians(1.623)
t0 = time.time()

# ---------------- f_K: Taylor (A, 4000) in D̄; obere Struktur nahe L (flint, wie v49)
a = np.array([complex(float(x[0] if isinstance(x, list) else x)) for x in
              json.loads((HIER.parent.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]])


def f_taylor(x):
    return np.polyval(a[::-1], x)


from flint import acb, arb, ctx  # noqa: E402
from kc_kneser import Kneser  # noqa: E402

KN = Kneser(200, laut=False)
R = KN.R
ctx.prec = 200
Lb = R.L
c = KN.c
pi = arb.pi()
TPI = 2 * pi * acb(0, 1)
tha = Lb.imag
sqL = (Lb / 2).exp()
codd = [c[n] if n % 2 == 1 else acb(0) for n in range(len(c))]


def mid(z):
    return acb(z.real.mid(), z.imag.mid())


def f_ober(x):
    x = acb(x.real, x.imag)
    w = mid(x - Lb)
    for _ in range(60):
        fv, df = R.P_dP_ball(w)
        dw = mid((fv - x) / df)
        w = mid(w - dw)
        if dw.abs_upper() < arb(2) ** (-180):
            break
    A = w.arg()
    while bool(A <= tha + 4 * pi):
        A += 2 * pi
    while bool(A > tha + 6 * pi):
        A -= 2 * pi
    zeta = acb(abs(w).log(), A) / Lb
    z = zeta - c[0]
    for _ in range(100):
        zn = mid(zeta - R._horner(c, (TPI * z).exp()))
        if (zn - z).abs_upper() < arb(2) ** (-180):
            z = zn
            break
        z = zn
    assert float(z.imag.mid()) > 0.36
    d = -2 * R._horner(codd, (TPI * z).exp())
    v = R.P_ball(sqL * w * (Lb * d).exp())
    return complex(float(v.real.mid()), float(v.imag.mid()))


def fK(x):
    """f_K am Punkt x (obere Halbebene oder reell); untere Hälfte über Konjugation."""
    if x.imag < 0:
        return fK(x.conjugate()).conjugate()
    if abs(x - L) < 0.35:
        return f_ober(x)
    if abs(x) <= rho * (1 + 1e-12):
        return complex(f_taylor(x))
    raise ValueError("ausserhalb")


# ---------------- 𝒞 (obere Hälfte, von ρ bis −ρ(1+s_L)), dann gespiegelt
def spirale(uA, r):
    return L + r * cmath.exp(1j * (math.radians(uA + 71.46) + TS * math.log(r / rho)))


def richtung(rad, seite):
    best = None
    for k in range(40001):
        phi = th - 0.2 + 0.4 * k / 40000
        if (seite == "linse" and phi <= th) or (seite == "gout" and phi >= th):
            continue
        t = rad * cmath.exp(1j * phi)
        dd = abs(abs(t - L) - R1)
        if best is None or dd < best[0]:
            best = (dd, phi)
    return best[1]


phi_g = richtung(rho, "gout")          # Γout-Punkt im Abstand r₁ von L
phi_l = richtung(rho * (1 + SL), "linse")
rr = np.geomspace(RMIN, R1, 160)
oben = []                               # Liste (Punkt, Art)
for phi in np.linspace(0, phi_g, 400):
    oben.append((rho * cmath.exp(1j * phi), "gout"))
g1 = rho * cmath.exp(1j * phi_g)
a_g = cmath.phase(g1 - L)
a2 = cmath.phase(spirale(UA2, R1) - L)
while a_g - a2 > math.pi:
    a_g -= 2 * math.pi
while a2 - a_g > math.pi:
    a_g += 2 * math.pi
for s in np.linspace(0, 1, 15)[1:-1]:
    oben.append((L + R1 * cmath.exp(1j * (a_g + (a2 - a_g) * s)), "verb"))
for r in rr[::-1]:
    oben.append((spirale(UA2, r), "sigma2"))
oben.append((L, "L"))
for r in rr:
    oben.append((spirale(UA1, r), "sigma1"))
a1 = cmath.phase(spirale(UA1, R1) - L)
l1 = rho * (1 + SL) * cmath.exp(1j * phi_l)
a_l = cmath.phase(l1 - L)
while a_l - a1 > math.pi:
    a_l -= 2 * math.pi
while a1 - a_l > math.pi:
    a_l += 2 * math.pi
for s in np.linspace(0, 1, 15)[1:-1]:
    oben.append((L + R1 * cmath.exp(1j * (a1 + (a_l - a1) * s)), "verb"))
for phi in np.linspace(phi_l, math.pi, 400):
    oben.append((rho * (1 + SL) * cmath.exp(1j * phi), "linse"))
pts = [p for p, _ in oben]
art = [k for _, k in oben]
# f_K auf der oberen Hälfte von 𝒞
fv = []
letzt = None
for p, k in oben:
    if k == "linse" and abs(p - L) >= 0.35:
        # f(t) = Log f(eᵗ) + 2πi m, m stetig
        e = cmath.exp(p)
        g = cmath.log(complex(f_taylor(e)) if abs(e) <= rho else fK(e))
        m = round(((letzt - g) / (2j * math.pi)).real) if letzt is not None else 0
        v = g + 2j * math.pi * m
    else:
        v = fK(p) if p != L else L
    fv.append(v)
    letzt = v
print(f"𝒞 obere Hälfte: {len(pts)} Punkte, f_K berechnet ({time.time() - t0:.0f}s)", flush=True)
# Kontrolle der Linse am Übergang (obere Struktur ↔ Log f(eᵗ)) bei |t − L| ≈ 0.3–0.4
for p, v, k in zip(pts, fv, art):
    if k == "linse" and 0.33 < abs(p - L) < 0.37:
        e = cmath.exp(p)
        g = cmath.log(complex(f_taylor(e)) if abs(e) <= rho else fK(e))
        print(f"  Übergangskontrolle t = {p:.4f}: obere Struktur {f_ober(p):.10f}, Log f(eᵗ) {g:.10f}", flush=True)
        break
C_pts = np.array(pts + [p.conjugate() for p in pts[::-1]])
C_f = np.array(fv + [v.conjugate() for v in fv[::-1]])
C_t = C_pts.copy()                       # Parameter t (für f(𝒞): f_K(f(t)) = eᵗ)
# D̄ ∖ U_{r₀}: Kreis plus Bögen
D_pts, D_f = [], []
for phi in np.linspace(-math.pi, math.pi, 1600, endpoint=False):
    z = rho * cmath.exp(1j * phi)
    if abs(z - L) >= R0 and abs(z - L.conjugate()) >= R0:
        D_pts.append(z)
for s in np.linspace(0, 1, 40):
    for Lx in (L, L.conjugate()):
        aa = cmath.phase(Lx) + math.pi / 2 + math.pi * s
        z = Lx + R0 * cmath.exp(1j * aa)
        if abs(z) <= rho:
            D_pts.append(z)
D_pts = np.array(D_pts)


def poly(points):
    return Polygon(np.c_[points.real, points.imag]).buffer(0)


circ = Polygon([(rho * math.cos(p), rho * math.sin(p)) for p in np.linspace(0, 2 * math.pi, 2000)])
P_D = circ.difference(shapely.Point(L.real, L.imag).buffer(R0, 256)).difference(shapely.Point(L.real, -L.imag).buffer(R0, 256))
P_C = poly(C_pts)
F_C = C_f                                 # Bildkurve f(𝒞)
P_F = poly(F_C)
E = shapely.union_all([P_D, P_C, P_F])
if E.geom_type != "Polygon":
    E = max(E.geoms, key=lambda g: g.area)
rand = E.exterior
print(f"E: Fläche {E.area:.3f}, Randlänge {rand.length:.3f}, Löcher {len(E.interiors)} ({time.time() - t0:.0f}s)", flush=True)
# Randpunkte mit f_K-Werten: Kandidaten aus allen drei Stücken, behalte die auf ∂E
kand_x, kand_f = [], []
for z, v in zip(C_pts, C_f):
    kand_x.append(z)
    kand_f.append(v)
for z, t in zip(F_C, C_t):
    kand_x.append(z)
    kand_f.append(cmath.exp(t))           # f_K(f(t)) = eᵗ
for z in D_pts:
    kand_x.append(z)
    kand_f.append(None)
kand_x = np.array(kand_x)
band = rand.buffer(1e-9)
auf = shapely.contains_xy(band, kand_x.real, kand_x.imag)
X, F = [], []
for z, v, ok in zip(kand_x, kand_f, auf):
    if ok:
        X.append(z)
        F.append(v if v is not None else fK(complex(z)))
X, F = np.array(X), np.array(F)
print(f"Randpunkte auf ∂E: {len(X)} (davon nahe L < 0.05: {int(np.sum(np.abs(X - L) < 0.05))}) ({time.time() - t0:.0f}s)", flush=True)
# Aussenwinkel bei L aus den Richtungen der Randpunkte nahe L
nah = X[(np.abs(X - L) > 1e-4) & (np.abs(X - L) < 2e-3)]
print("  Richtungen der Randpunkte bei |x − L| ∈ [1e-4, 2e-3]: " +
      f"{np.degrees(np.angle(nah - L)).min():.1f}° … {np.degrees(np.angle(nah - L)).max():.1f}° (Anzahl {len(nah)})", flush=True)
alpha_E = (UA1 - UA2 - 38.3) / 180
print(f"Vorhersage: Aussenwinkel {UA1 - UA2 - 38.3:.1f}° ⇒ s_E = 5.447·α_E/(1+T²) = {5.447 * alpha_E / (1 + 0.02833 ** 2):.3f}", flush=True)
# ---------------- diskrete Bestapproximation (Arnoldi + Lawson)
c0 = complex(E.centroid.x, E.centroid.y)
Rs = np.max(np.abs(X - c0))
Z = (X - c0) / Rs
erg = {}
for N in (20, 40, 60, 80, 120, 160, 200, 250, 300):
    m = len(Z)
    Q = np.zeros((m, N + 1), dtype=complex)
    Q[:, 0] = 1
    for k in range(1, N + 1):
        q = Z * Q[:, k - 1]
        for _ in range(2):                                   # zweimal Gram–Schmidt
            h = (Q[:, :k].conj().T @ q) / m
            q = q - Q[:, :k] @ h
        Q[:, k] = q / (np.linalg.norm(q) / math.sqrt(m))
    wgt = np.ones(m) / m
    best = None
    for it in range(40):
        sw = np.sqrt(wgt)
        coef, *_ = np.linalg.lstsq(Q * sw[:, None], F * sw, rcond=None)
        err = np.abs(Q @ coef - F)
        if best is None or err.max() < best[0]:
            best = (err.max(), np.argmax(err))
        wgt = wgt * err
        wgt /= wgt.sum()
    erg[N] = best[0]
    zmax = X[best[1]]
    print(f"N = {N}: max Fehler ≈ {best[0]:.3e} (bei x = {zmax:.4f}, |x − L| = {abs(zmax - L):.2e}) ({time.time() - t0:.0f}s)", flush=True)
Ns = sorted(erg)
for i in range(1, len(Ns)):
    s_loc = -math.log(erg[Ns[i]] / erg[Ns[i - 1]]) / math.log(Ns[i] / Ns[i - 1])
    print(f"  lokaler Exponent {Ns[i - 1]}→{Ns[i]}: {s_loc:.2f}")
json.dump(dict(uA1=UA1, uA2=UA2, r0=R0, fehler={str(k): v for k, v in erg.items()}), open(HIER / f"v52_sE_{int(UA1)}_{int(UA2)}_r1{R1}.json", "w", encoding="utf-8"))
