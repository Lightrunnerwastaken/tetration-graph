"""Strang V — F3 (HAUPTSATZ_KR0 Fassung 2 §2.4; K BRETT 11:05/11:16): sup_𝒞 |p_N − f_K| am echten p_N bis N = 2000.

𝒞 wie in v52/v51: Σ₁ (u_A = 80°), Σ₂ (u_A = −90°) bis r₁ = 0.05, Verbindungsbögen |t − L| = r₁, Linse |t| = ρ(1 + 0.03) über Γin, Γout-Bogen |t| = ρ (obere Hälfte;
unten konjugiert, p_N und f_K reell). f_K: Taylor (A, 4000) in D̄, Ps obere Struktur nahe L (|t − L| < 0.35), auf der Linse fern von L f(t) = Log f(eᵗ) (stetiger Zweig).
p_N: A (N ≤ 1000, 60 Stellen), v33 (1500, 2000); Horner am Mittelpunkt mit 600 Bit.
Ausgabe je N: max|p_N − f_K| je Teil (Linse, Σ₁, Σ₂, Verbindung, Γout) mit Lage, dazu getrennt die Linse am Steg (180° ± 30°).
Aufruf: python v53_F3.py [N-Liste]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "pv_kopie"))
sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "200,400,600,800,1000,1500,2000").split(",")]
L = complex(0.31813150520476413, 1.3372357014306895)
rho, th = abs(L), cmath.phase(L)
TS = math.radians(1.623)
R1, SL, RMIN = 0.05, 0.03, 1e-5
t0 = time.time()
a = np.array([complex(float(x[0] if isinstance(x, list) else x)) for x in
              json.loads((HIER.parent.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]])

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
    d = -2 * R._horner(codd, (TPI * z).exp())
    v = R.P_ball(sqL * w * (Lb * d).exp())
    return complex(float(v.real.mid()), float(v.imag.mid()))


def spirale(uA, r):
    return L + r * cmath.exp(1j * (math.radians(uA + 71.46) + TS * math.log(r / rho)))


def richtung(rad, seite):
    best = None
    for k in range(40001):
        phi = th - 0.2 + 0.4 * k / 40000
        if (seite == "linse" and phi <= th) or (seite == "gout" and phi >= th):
            continue
        dd = abs(abs(rad * cmath.exp(1j * phi) - L) - R1)
        if best is None or dd < best[0]:
            best = (dd, phi)
    return best[1]


phi_g, phi_l = richtung(rho, "gout"), richtung(rho * (1 + SL), "linse")
rr = np.geomspace(RMIN, R1, 60)
pkt = []
for phi in np.linspace(0, phi_g, 150):
    pkt.append((rho * cmath.exp(1j * phi), "gout", math.degrees(phi)))
for r in rr:
    pkt.append((spirale(-90, r), "sigma2", r))
    pkt.append((spirale(80, r), "sigma1", r))
a1, a2 = cmath.phase(spirale(80, R1) - L), cmath.phase(spirale(-90, R1) - L)
a_l, a_g = cmath.phase(rho * (1 + SL) * cmath.exp(1j * phi_l) - L), cmath.phase(rho * cmath.exp(1j * phi_g) - L)
for s in np.linspace(0, 1, 12):
    pkt.append((L + R1 * cmath.exp(1j * (a1 + (a_l - a1) * s)), "verb", None))
    pkt.append((L + R1 * cmath.exp(1j * (a2 + (a_g - a2) * s)), "verb", None))
linse = []
for phi in np.linspace(phi_l, math.pi, 300):
    linse.append(rho * (1 + SL) * cmath.exp(1j * phi))
# f_K-Werte
werte = []
for p, k, info in pkt:
    werte.append(f_ober(p) if abs(p - L) < 0.35 else complex(np.polyval(a[::-1], p)))
letzt = None
for p in linse:
    if abs(p - L) < 0.35:
        v = f_ober(p)
    else:
        e = cmath.exp(p)
        g = cmath.log(complex(np.polyval(a[::-1], e)))
        m = round(((letzt - g) / (2j * math.pi)).real)
        v = g + 2j * math.pi * m
    pkt.append((p, "linse", math.degrees(cmath.phase(p))))
    werte.append(v)
    letzt = v
print(f"𝒞: {len(pkt)} Punkte (obere Hälfte), f_K berechnet ({time.time() - t0:.0f}s)", flush=True)


def lade(N):
    if N <= 1000:
        d = json.loads((HIER.parent.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"][str(N)]
    else:
        d = json.loads((HIER / "sg" / f"v33_pN_{N}.json").read_text(encoding="utf-8"))["p"]
    ctx.prec = 600
    cf = [arb(s) for s in d]
    return cf


erg = {}
for N in NS:
    try:
        cf = lade(N)
    except (FileNotFoundError, KeyError):
        print(f"N = {N}: keine Daten")
        continue
    ctx.prec = 600
    teile = {}
    for (p, k, info), v in zip(pkt, werte):
        x = acb(p.real, p.imag)
        s = acb(0)
        for co in reversed(cf):
            s = s * x + co
        e = abs(complex(float(s.real.mid()), float(s.imag.mid())) - v)
        key = k
        if k == "linse" and info is not None and 150 <= info <= 180:
            key = "steg"
        if key not in teile or e > teile[key][0]:
            teile[key] = (e, info if info is not None else abs(p - L))
    ctx.prec = 200
    erg[N] = {k: v[0] for k, v in teile.items()}
    print(f"N = {N} ({time.time() - t0:.0f}s): " + ", ".join(f"{k} {v[0]:.2e} (bei {v[1]:.4g})" for k, v in sorted(teile.items())), flush=True)
json.dump({str(k): v for k, v in erg.items()}, open(HIER / "v53_F3.json", "w", encoding="utf-8"), indent=0)
