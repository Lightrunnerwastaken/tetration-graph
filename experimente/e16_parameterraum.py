"""exp.parameterraum — Regionen der Basen-Ebene (Bild 1): Periode des Turms b^{b^{…}} über komplexen b.

 Parameter t = −W_0(−Log b) (= Multiplikator des Hauptfixpunkts), b = exp(t e^{−t}).
 (1) Konforme Parametrisierung: t ↦ b(t) ist auf |t| < 1 umkehrbar durch −W_0(−Log b).
 (2) Periode 1 ⇔ |t| < 1 (beide Richtungen, Gitter in der b-Ebene).
 (3) Reelle Achse: (0, e^{−e}) Periode 2, (e^{−e}, η) Periode 1, b > η Flucht.
 (4) Basis i (Justus, Complex Tetration Dynamics, 2026): L, λ, |λ|, arg λ; Spirale: (z_{n+1} − L)/(z_n − L) → λ.
 (5) Knospen: der q-Zyklus umläuft den (abstossend gewordenen) Fixpunkt mit Drehzahl p/q.
 (6) Rand, rationaler Winkel: parabolisch, |z_n − L| ≈ C n^{−1/q} (q = 2: b = e^{−e}; q = 3).
"""
import cmath
import math

import numpy as np
from scipy.special import lambertw

from gemeinsam import Lauf, plt

lauf = Lauf("exp.parameterraum", __file__)
ETA, EE = math.exp(1 / math.e), math.exp(-math.e)


def t_von_b(b):
    return -lambertw(-np.log(b), 0)


def b_von_t(t):
    return np.exp(t * np.exp(-t))


def periode(b, n_vor=3000, qmax=12, tol=1e-8):
    u = np.log(b)
    z = np.ones_like(b)
    with np.errstate(all="ignore"):
        for _ in range(n_vor):
            z = np.exp(u * z)
            z = np.where(np.isfinite(z) & (np.abs(z) < 1e12), z, np.nan)
        bahn = [z]
        for _ in range(qmax):
            bahn.append(np.exp(u * bahn[-1]))
    per = np.zeros(b.shape, dtype=int)
    for q in range(qmax, 0, -1):
        per = np.where(np.abs(bahn[q] - bahn[0]) < tol * (1 + np.abs(bahn[0])), q, per)
    return np.where(np.isnan(bahn[0]), -1, per)        # −1 Flucht, 0 unklar (kein Zyklus ≤ qmax erkannt)


# (1) Parametrisierung
rng = np.random.default_rng(7)
tt = 0.95 * np.sqrt(rng.random(4000)) * np.exp(2j * np.pi * rng.random(4000))
fehler = np.max(np.abs(t_von_b(b_von_t(tt)) - tt))
lauf.pruefe("t ↦ exp(t e^{−t}) auf |t| < 0.95 umkehrbar durch −W_0(−Log b) (4000 Zufallspunkte)", fehler < 1e-12, f"{fehler:.1e}")

# (2) Gitter in der b-Ebene
re = np.linspace(-1.0, 3.0, 520)
im = np.linspace(-2.4, 2.4, 520)
B = re[None, :] + 1j * im[:, None]
B[B == 0] = 1e-9
P = periode(B)
T = np.abs(t_von_b(B))
innen = T < 0.97
lauf.pruefe(f"Gegenbeispielsuche: alle {int(innen.sum())} Gitterbasen mit |t| < 0.97 haben Periode 1", int(np.sum(innen & (P != 1))) == 0,
            f"{int(np.sum(innen & (P != 1)))} Gegenbeispiele")
aussen1 = (P == 1) & (T > 1.02)
lauf.pruefe("Umkehrung: keine Gitterbasis mit |t| > 1.02 konvergiert (anziehender Fixpunkt nur auf dem Hauptzweig)", int(aussen1.sum()) == 0,
            f"{int(aussen1.sum())} von {int((T > 1.02).sum())}")
anteile = {str(k): int(np.sum(P == k)) for k in range(-1, 13)}

# (3) reelle Achse
reell = {"(0, e^-e)": np.linspace(0.005, EE - 0.004, 40), "(e^-e, eta)": np.linspace(EE + 0.004, ETA - 0.004, 60),
         "(eta, 4)": np.linspace(ETA + 0.01, 4, 40)}
pr = {k: periode(v.astype(complex), n_vor=20000) for k, v in reell.items()}
ok = np.all(pr["(0, e^-e)"] == 2) and np.all(pr["(e^-e, eta)"] == 1) and np.all(pr["(eta, 4)"] == -1)
lauf.pruefe("reelle Achse: Periode 2 auf (0, e^{−e}), 1 auf (e^{−e}, η), Flucht über η", bool(ok),
            {k: sorted(set(int(x) for x in v)) for k, v in pr.items()})

# (4) Basis i
u = 1j * math.pi / 2
L = -complex(lambertw(-u, 0)) / u
lam = L * u
notiz = {"L": 0.438282936727 + 0.360592471871j, "lam": -0.566417330285 + 0.688453227108j, "betrag": 0.891513566, "arg": 2.259249554}
abw = max(abs(L - notiz["L"]), abs(lam - notiz["lam"]), abs(abs(lam) - notiz["betrag"]), abs(cmath.phase(lam) - notiz["arg"]))
lauf.pruefe("Basis i: Fixpunkt, Multiplikator, |λ|, arg λ wie in der Notiz (12 Stellen)", abw < 5e-10, f"{abw:.1e}")
z, quot = 1 + 0j, []
for n in range(400):
    z_neu = cmath.exp(u * z)
    if n >= 20:
        quot.append((z_neu - L) / (z - L))
    z = z_neu
spirale = abs(quot[150] - lam)
lauf.pruefe("Spirale: (z_{n+1} − L)/(z_n − L) → λ (Schrumpfung |λ|, Drehung arg λ je Schritt)", spirale < 1e-6, f"{spirale:.1e} bei n = 170")

# (5) Knospen-Drehzahl
knospen, treffer = [], 0
for q in range(2, 8):
    for p in range(1, q):
        if math.gcd(p, q) != 1:
            continue
        t = (1 + 0.5 / q ** 2) * cmath.exp(2j * math.pi * p / q)
        uu = t * cmath.exp(-t)
        Lq = cmath.exp(t)
        z = 1 + 0j
        for _ in range(60000):
            z = cmath.exp(uu * z)
        zyk = [z]
        for _ in range(q - 1):
            zyk.append(cmath.exp(uu * zyk[-1]))
        geschlossen = abs(cmath.exp(uu * zyk[-1]) - zyk[0]) < 1e-8
        ordnung = sorted(range(q), key=lambda j: cmath.phase(zyk[j] - Lq) % (2 * math.pi))
        rang = {j: k for k, j in enumerate(ordnung)}
        schritte = {(rang[(j + 1) % q] - rang[j]) % q for j in range(q)}
        drehzahl = schritte.pop() if len(schritte) == 1 else None
        knospen.append({"p": p, "q": q, "geschlossen": geschlossen, "drehzahl": drehzahl})
        treffer += geschlossen and drehzahl == p
lauf.pruefe("Knospen p/q (q ≤ 7): q-Zyklus umläuft den Fixpunkt mit Drehzahl p/q", treffer == len(knospen), f"{treffer}/{len(knospen)}")

# (6) parabolische Randpunkte
exponenten = {}
for q in (2, 3):
    t = cmath.exp(2j * math.pi / q)
    uu = t * cmath.exp(-t)
    Lq = cmath.exp(t)
    z, mess = 1 + 0j, {}
    for n in range(1, 400001):
        z = cmath.exp(uu * z)
        if n in (20000, 400000):
            mess[n] = abs(z - Lq)
    exponenten[q] = -math.log(mess[400000] / mess[20000]) / math.log(20)
lauf.pruefe("parabolischer Rand: |z_n − L| ≈ C n^{−1/q} für q = 2 (b = e^{−e}) und q = 3 (Exponent ± 0.05)",
            all(abs(exponenten[q] - 1 / q) < 0.05 for q in exponenten), {q: round(v, 4) for q, v in exponenten.items()})
lauf.daten = {"perioden_anteile": anteile, "knospen": knospen, "exponenten": exponenten, "basis_i": {"L": str(L), "lam": str(lam)},
              "gitter": [len(re), len(im)]}

# Bild
fig, ax = plt.subplots(1, 2, figsize=(11, 5.2), gridspec_kw={"width_ratios": [1.25, 1]})
farbe = np.where(P > 0, P, np.nan).astype(float)
ax[0].imshow(farbe, extent=[re[0], re[-1], im[0], im[-1]], origin="lower", cmap=plt.get_cmap("viridis", 12), vmin=0.5, vmax=12.5,
             interpolation="nearest")
ax[0].imshow(np.where(P == -1, 1.0, np.nan), extent=[re[0], re[-1], im[0], im[-1]], origin="lower", cmap="Greys", vmin=0, vmax=3,
             interpolation="nearest")
th = np.linspace(0, 2 * np.pi, 800)
rand = b_von_t(np.exp(1j * th))
ax[0].plot(rand.real, rand.imag, color="#b5523b", lw=1.1, label="Shell–Thron-Rand b(θ) = exp(e^{iθ − e^{iθ}})")
for p, q in ((1, 2), (1, 3), (1, 4), (1, 5), (2, 5)):
    bk = b_von_t(1.08 * np.exp(2j * np.pi * p / q))
    ax[0].annotate(f"{p}/{q}", (bk.real, bk.imag), fontsize=7, color="#222")
ax[0].plot([EE, ETA], [0, 0], "o", ms=3, color="#b5523b")
ax[0].set_xlabel("Re b"); ax[0].set_ylabel("Im b")
ax[0].set_title("Periode des Turms (Farbe), grau = Flucht, weiss = unklar")
ax[0].legend(fontsize=7, loc="lower left")
z, orbit = 1 + 0j, [1 + 0j]
for _ in range(60):
    z = cmath.exp(u * z)
    orbit.append(z)
O = np.array(orbit)
ax[1].plot(O.real, O.imag, ".-", lw=0.6, ms=4, color="#3d6fb6", label="z_{n+1} = i^{z_n}, z_0 = 1")
ax[1].plot([L.real], [L.imag], "o", color="#b5523b", label=f"L, λ = {abs(lam):.4f}·e^({cmath.phase(lam):.4f} i)")
ax[1].set_aspect("equal"); ax[1].legend(fontsize=8); ax[1].set_title("Basis i: Spirale mit Drehung arg λ und Schrumpfung |λ|")
lauf.bild(fig)
lauf.zusammenfassung = (f"Periode 1 genau für |t| < 1 (0 Gegenbeispiele in beide Richtungen auf {len(re)}×{len(im)} Gitter); reelle Achse 2/1/Flucht; "
                        f"Basis i wie in der Notiz; Knospen-Drehzahl p/q in {treffer}/{len(knospen)}; parabolische Exponenten "
                        + ", ".join(f"q={q}: {v:.3f}" for q, v in exponenten.items()) + ".")
lauf.schreibe()
