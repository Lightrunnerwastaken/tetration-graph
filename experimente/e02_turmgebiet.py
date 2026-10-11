"""exp.turmgebiet — Konvergenz des unendlichen Turms über komplexen Basen.

Parameter c = ln b. Iteration z ↦ exp(c z) ab z = 1. Innerhalb des Shell–Thron-Gebiets
(c = t e^{−t}, |t| < 1) muss der Turm konvergieren (Gegenbeispielsuche), in den Knospen
bei t = r e^{2πi p/q}, r knapp > 1, sollte ein q-Zyklus entstehen (Vermutung fp.verm.knospen).
"""
import cmath
import math

import numpy as np
from scipy.special import lambertw

from gemeinsam import Lauf, plt

lauf = Lauf("exp.turmgebiet", __file__)


def periode(c, n_vor=4000, qmax=12, tol=1e-9):
    z = np.ones_like(c)
    with np.errstate(all="ignore"):
        for _ in range(n_vor):
            z = np.exp(c * z)
            z = np.where(np.abs(z) > 1e12, np.nan, z)
        bahn = [z]
        for _ in range(qmax):
            bahn.append(np.exp(c * bahn[-1]))
    per = np.zeros(c.shape, dtype=int)
    for q in range(qmax, 0, -1):
        gleich = np.abs(bahn[q] - bahn[0]) < tol * (1 + np.abs(bahn[0]))
        per = np.where(gleich, q, per)
    per = np.where(np.isnan(bahn[0]), -1, per)
    return per


# Gitter in der c-Ebene
re = np.linspace(-3.2, 0.6, 420)
im = np.linspace(-2.2, 2.2, 420)
C = re[None, :] + 1j * im[:, None]
P = periode(C)
t = -lambertw(-C, 0)
innen = np.abs(t) < 0.97
gegen = int(np.sum(innen & (P != 1)))
lauf.pruefe(f"Gegenbeispielsuche: alle {int(innen.sum())} Gitterpunkte mit |t| < 0.97 konvergieren", gegen == 0,
            f"{gegen} Gegenbeispiele")

# Knospen
knospen = []
for q in range(2, 8):
    for p in range(1, q):
        if math.gcd(p, q) != 1:
            continue
        for r in (1 + 0.3 / q ** 2, 1 + 0.8 / q ** 2):
            tt = r * cmath.exp(2j * math.pi * p / q)
            c = tt * cmath.exp(-tt)
            per = int(periode(np.array([c]))[0])
            knospen.append({"p": p, "q": q, "r": round(r, 5), "periode": per})
treffer = [k for k in knospen if k["periode"] == k["q"]]
lauf.pruefe("Knospen: Periode = q in allen getesteten Fällen", len(treffer) == len(knospen),
            f"{len(treffer)}/{len(knospen)}")
lauf.daten = {"knospen": knospen, "gitter": [len(re), len(im)], "gegenbeispiele_innen": gegen}

fig, ax = plt.subplots(figsize=(6.5, 5.5))
farbe = np.where(P < 0, 0, P).astype(float)
farbe[P == 0] = np.nan
cm = plt.get_cmap("viridis", 12)
ax.imshow(np.where(P > 0, P, np.nan), extent=[re[0], re[-1], im[0], im[-1]], origin="lower", cmap=cm, vmin=0.5, vmax=12.5,
          interpolation="nearest")
phi = np.linspace(0, 2 * np.pi, 600)
tt = np.exp(1j * phi)
rand = tt * np.exp(-tt)
ax.plot(rand.real, rand.imag, color="#b5523b", lw=1.2, label="Shell–Thron-Rand  c = t e^{−t}, |t| = 1")
ax.set_xlabel("Re ln b")
ax.set_ylabel("Im ln b")
ax.set_title("Grenzperiode des Turms b^{b^{…}} (Farbe = Periode, weiss = Flucht/unklar)")
ax.legend(fontsize=8, loc="lower left")
lauf.bild(fig)
lauf.zusammenfassung = (f"{int(innen.sum())} Basen mit |t|<0.97: alle konvergent (0 Gegenbeispiele). "
                        f"Knospen p/q, q ≤ 7: Periode q in {len(treffer)}/{len(knospen)} Fällen.")
lauf.schreibe()
