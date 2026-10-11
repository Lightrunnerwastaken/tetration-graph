"""exp.mittelwerte — Mittelwerte über Basen auf und in der Shell–Thron-Region (Justus, Complex Tetration Dynamics, 2026).

Basen b(w) = exp(w e^{−w}) (Log b fest), Fixpunkt e^w, Multiplikator w.
 (1) Satz 1 der Notiz: (1/2π) ∫ T_{b(e^{iθ})}(n) dθ = 1 für ganze n; auch für |w| = r < 1.
 (2) Korollar 2: T_{z^{1/z}}(n) = (1/2π) ∫ T_{b(θ)}(n) / (1 − ln z · e^{−iθ}) dθ, geprüft bei z = 2 (Basis √2).
 (3) Vermutung der Notiz für nichtganze Höhen: (1/2π) ∫ T_{b(re^{iθ})}(t) dθ = 1 auch für t = 1/2?
     Gegenbeispielsuche mit regulärer (Koenigs-)Tetration, zwei Zweige von λ^t.
 (4) Siegel-Fall: Cesàro-Mittel des Orbits bei goldenem Rotationswinkel → Fixpunkt.
"""
import cmath
import math

import mpmath as mp
import numpy as np

from gemeinsam import Lauf, plt
from koenigs_komplex import KoenigsKomplex

mp.mp.dps = 30
lauf = Lauf("exp.mittelwerte", __file__)


def Z(n, w):
    u = w * mp.exp(-w)
    z = mp.mpc(1)
    for _ in range(n):
        z = mp.exp(z * u)
    return z


def kreismittel(f, r, M=256):
    return mp.fsum(f(r * mp.expjpi(2 * mp.mpf(j) / M)) for j in range(M)) / M


# (1) ganze Höhen
for r in (mp.mpf(1), mp.mpf("0.6")):
    abw = max(abs(kreismittel(lambda w: Z(n, w), r, M=2048 if r == 1 else 256) - 1) for n in (1, 2, 3, 5, 8))
    lauf.pruefe(f"Satz 1: Mittel von T_b(w)(n) über |w| = {mp.nstr(r, 2)} ist 1 (n = 1, 2, 3, 5, 8)", abw < mp.mpf(10) ** -25, mp.nstr(abw, 3))

# (2) Korollar 2 bei z = 2
z = mp.mpf(2)
for n in (2, 3):
    links = Z(n, mp.log(z))                                   # Basis exp(ln z · e^{−ln z}) = z^{1/z} = √2
    rechts = mp.fsum(Z(n, mp.expjpi(2 * mp.mpf(j) / 256)) / (1 - mp.log(z) * mp.expjpi(-2 * mp.mpf(j) / 256)) for j in range(256)) / 256
    lauf.pruefe(f"Korollar 2 bei z = 2, n = {n}: T_√2(n) aus dem Randintegral", abs(links - rechts) < mp.mpf(10) ** -20,
                f"{mp.nstr(mp.re(links), 12)} vs {mp.nstr(rechts, 12)}")

# (3) nichtganze Höhe t = 1/2 im Inneren (|w| = r < 1)
M = 48
ergebnis_bruch = {}
for r in (mp.mpf("0.5"), mp.mpf("0.8")):
    for name, zweig in (("Hauptzweig (−π, π]", lambda th: th if th <= mp.pi else th - 2 * mp.pi), ("Zweig [0, 2π)", lambda th: th)):
        werte = []
        for j in range(M):
            th = 2 * mp.pi * j / M
            w = r * mp.expj(th)
            K = KoenigsKomplex(w, N=30, log_lam=mp.mpc(mp.log(r), zweig(th)))
            werte.append(K.tet(mp.mpf("0.5"), schritte=12))
        mittel = mp.fsum(werte) / M
        ergebnis_bruch[f"r={mp.nstr(r, 2)}, {name}"] = mp.nstr(mittel, 10)
        print("  t=1/2", r, name, mp.nstr(mittel, 12))
abweichung = max(abs(mp.mpc(v) - 1) for v in map(mp.mpmathify, ergebnis_bruch.values()))
lauf.pruefe("Gegenbeispiel: für t = 1/2 ist das Mittel nicht 1 (mindestens ein Fall weicht > 1e-3 ab)", abweichung > 1e-3,
            mp.nstr(abweichung, 4))
zweige = list(ergebnis_bruch.values())
lauf.pruefe("das Mittel hängt vom Zweig von λ^t ab (r = 0.5)", abs(mp.mpmathify(zweige[0]) - mp.mpmathify(zweige[1])) > 1e-4,
            f"{zweige[0]} vs {zweige[1]}")

# (4) Siegel-Fall: goldener Rotationswinkel
phi = (math.sqrt(5) - 1) / 2
w = cmath.exp(2j * math.pi * phi)
u = w * cmath.exp(-w)
L = cmath.exp(w)
zz = 1 + 0j
summe, bahn, fehler = 0j, [], {}
for n in range(1, 200001):
    zz = cmath.exp(zz * u)
    summe += zz
    if n % 4 == 0 and n <= 40000:
        bahn.append(zz)
    if n in (1000, 10000, 100000, 200000):
        fehler[n] = abs(summe / n - L)
lauf.pruefe("Siegel-Fall (goldener Winkel): Cesàro-Mittel → Fixpunkt, Fehler fällt etwa wie 1/N",
            fehler[200000] < 1e-3 and fehler[200000] < fehler[1000] / 30, {k: f"{v:.2e}" for k, v in fehler.items()})
radien = np.abs(np.array(bahn) - L)
lauf.daten = {"bruch_t_halbe": ergebnis_bruch, "siegel_fehler": {k: v for k, v in fehler.items()},
              "siegel_radius_min_max": [float(radien.min()), float(radien.max())]}

fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
B = np.array(bahn)
ax[0].plot(B.real, B.imag, ",", color="#3d6fb6", alpha=0.6)
ax[0].plot([L.real], [L.imag], "o", color="#b5523b", label="Fixpunkt e^w")
ax[0].plot([(summe / 200000).real], [(summe / 200000).imag], "x", color="#2f8f83", ms=10, label="Cesàro-Mittel (N = 2·10⁵)")
ax[0].set_aspect("equal"); ax[0].legend(fontsize=8)
ax[0].set_title("Rotationsbasis, goldener Winkel: Orbit auf invarianter Kurve")
ns = sorted(fehler)
ax[1].loglog(ns, [fehler[n] for n in ns], "o-", label="|Mittel − Fixpunkt|")
ax[1].loglog(ns, [fehler[ns[0]] * ns[0] / n for n in ns], "--", color="#888", label="∝ 1/N")
ax[1].set_xlabel("N"); ax[1].legend(fontsize=8); ax[1].set_title("Mittelwert konvergiert gegen den Fixpunkt")
lauf.bild(fig)
lauf.zusammenfassung = ("Satz 1 (ganze Höhen) und Korollar 2 bestätigt. Für t = 1/2 weicht das Kreismittel von 1 ab und hängt vom Zweig von λ^t ab: "
                        + "; ".join(f"{k}: {v}" for k, v in ergebnis_bruch.items())
                        + f". Siegel-Fall: Cesàro-Fehler {fehler[200000]:.1e} bei N = 2·10⁵.")
lauf.schreibe()
