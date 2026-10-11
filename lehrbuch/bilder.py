"""Abbildungen für das Lehrbuch (Kapitel 0 und Grundlagen). Ausgabe: lehrbuch/bilder/<name>.svg

    python lehrbuch/bilder.py
Kneser-Werte kommen aus dem Zwischenspeicher der Experimente (gemeinsam.kneser).
"""
import math
import sys
from pathlib import Path

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "experimente"))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

from gemeinsam import kneser, plt  # noqa: E402
from regulaer import Regulaer  # noqa: E402

plt.rcParams.update({"svg.fonttype": "none", "font.size": 10})
AUS = HIER / "bilder"
AUS.mkdir(exist_ok=True)
E = math.e
ETA = E ** (1 / E)
ROT, BLAU, GRUEN, VIOLETT, GRAU = "#b5523b", "#3d6fb6", "#2f8f83", "#7a49c2", "#8a93a3"


def speichern(fig, name):
    fig.tight_layout()
    fig.savefig(AUS / f"{name}.svg", format="svg", metadata={"Date": None})
    plt.close(fig)
    print("  ", name)


# ------------------------------------------------------------------ k0_rechenarten: Addition, Multiplikation, Potenz, Turm
fig, ax = plt.subplots(figsize=(6.4, 3.6))
n = np.arange(1, 5)
reihen = {"2 + 2 + … (n Summanden) = 2n": 2 * n, "2 · 2 · … (n Faktoren) = 2ⁿ": 2.0 ** n,
          "Turm 2^2^…^2 (n Zweien) = 2↑↑n": np.array([2, 4, 16, 65536], dtype=float)}
for (label, y), c, m in zip(reihen.items(), (GRAU, BLAU, ROT), ("o", "s", "D")):
    ax.semilogy(n, y, marker=m, color=c, label=label)
for x, y in zip(n, [2, 4, 16, 65536]):
    ax.annotate(f"{y:,}".replace(",", " "), (x, y), textcoords="offset points", xytext=(6, -2), fontsize=9, color=ROT)
ax.text(3.05, 2.5, "2↑↑5 hätte schon 19 729 Stellen", fontsize=9, color=ROT)
ax.set_xticks(n); ax.set_xlabel("n"); ax.set_ylabel("Wert (logarithmische Achse)")
ax.set_title("Vier Stufen derselben Idee: wiederholen")
ax.legend(fontsize=8, loc="upper left")
speichern(fig, "k0_rechenarten")

# ------------------------------------------------------------------ k0_treppe: b^x gegen x, Treppe zum Fixpunkt
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4))
xs = np.linspace(-0.5, 5, 400)
ax[0].plot(xs, xs, color="#222", lw=1, label="y = x")
for b, c in ((1.3, GRUEN), (ETA, BLAU), (1.6, ROT)):
    ax[0].plot(xs, b ** xs, color=c, label=f"y = {b:.4g}^x" if b != ETA else "y = η^x  (η ≈ 1.4447)")
ax[0].plot([E], [E], "o", color=BLAU, ms=5)
ax[0].annotate("berührt bei x = e", (E, E), xytext=(3.0, 1.2), fontsize=8, arrowprops=dict(arrowstyle="-", color=GRAU))
ax[0].set_xlim(-0.5, 5); ax[0].set_ylim(-0.5, 6); ax[0].set_aspect("equal")
ax[0].set_title("Wo schneidet y = bˣ die Gerade y = x?"); ax[0].legend(fontsize=8, loc="upper left")
b = math.sqrt(2)
xs = np.linspace(0, 4.5, 300)
ax[1].plot(xs, xs, color="#222", lw=1)
ax[1].plot(xs, b ** xs, color=BLAU, label="y = (√2)^x")
z = 1.0
for _ in range(14):
    z2 = b ** z
    ax[1].plot([z, z, z2], [z, z2, z2], color=ROT, lw=0.9)
    z = z2
ax[1].plot([2, 4], [2, 4], "o", color="#222", ms=4)
ax[1].annotate("a = 2 (zieht an)", (2, 2), xytext=(2.3, 1.2), fontsize=8)
ax[1].annotate("A = 4 (stösst ab)", (4, 4), xytext=(2.6, 4.25), fontsize=8)
ax[1].set_xlim(0, 4.5); ax[1].set_ylim(0, 4.5); ax[1].set_aspect("equal")
ax[1].set_title("Turm 1, √2, √2^√2, … : eine Treppe zu 2"); ax[1].legend(fontsize=8, loc="upper left")
speichern(fig, "k0_treppe")

# ------------------------------------------------------------------ k0_kleine_basen: Spirale (b = 0.3) und 2-Zyklus (b = 0.03)
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4))
for a_, b, titel in ((ax[0], 0.3, "b = 0.3: hin und her, aber zur Ruhe"), (ax[1], 0.03, "b = 0.03: ewiges Springen (2-Zyklus)")):
    xs = np.linspace(0, 1.05, 300)
    a_.plot(xs, xs, color="#222", lw=1)
    a_.plot(xs, b ** xs, color=BLAU, label=f"y = {b}^x")
    z = 1.0
    for _ in range(30):
        z2 = b ** z
        a_.plot([z, z, z2], [z, z2, z2], color=ROT, lw=0.8)
        z = z2
    a_.set_xlim(0, 1.05); a_.set_ylim(0, 1.05); a_.set_aspect("equal"); a_.set_title(titel); a_.legend(fontsize=8)
speichern(fig, "k0_kleine_basen")

# ------------------------------------------------------------------ k0_xhochx: x^x und die Türme von x
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4))
xs = np.linspace(1e-4, 2, 500)
ax[0].plot(xs, xs ** xs, color=BLAU, label="y = xˣ")
ax[0].plot([1 / E], [(1 / E) ** (1 / E)], "o", color=ROT)
ax[0].annotate("Minimum bei x = 1/e,\nWert e^(−1/e) ≈ 0.692", (1 / E, (1 / E) ** (1 / E)), xytext=(0.55, 1.35), fontsize=8,
               arrowprops=dict(arrowstyle="-", color=GRAU))
ax[0].set_ylim(0, 4.2); ax[0].set_xlabel("x"); ax[0].set_title("Die Selbstpotenz xˣ"); ax[0].legend(fontsize=8)
xs = np.linspace(1e-4, ETA, 600)
t = xs.copy()
for n in range(1, 9):
    if n > 1:
        t = xs ** t
    ax[1].plot(xs, t, lw=0.9, color=plt.cm.viridis(n / 9), label=f"{n}x" if n in (1, 2, 3, 8) else None)
ax[1].axvline(math.exp(-E), ls=":", color=GRAU)
ax[1].axvline(ETA, ls=":", color=GRAU)
ax[1].text(math.exp(-E) + 0.01, 2.6, "e^(−e)", fontsize=8, color=GRAU)
ax[1].text(ETA - 0.2, 2.6, "η", fontsize=8, color=GRAU)
ax[1].set_ylim(0, 2.9); ax[1].set_xlabel("x"); ax[1].set_title("Türme ⁿx = x^x^…^x (n Stück)")
ax[1].legend(fontsize=8, title="Turmhöhe n", title_fontsize=8, loc="upper left")
speichern(fig, "k0_xhochx")

# ------------------------------------------------------------------ k0_halbe_hoehe: Kneser T_e, lineare Tetration, reguläre √2-Tetration
mp.mp.dps = 30
xs = [mp.mpf(i) / 20 for i in range(-38, 47)]
grund = [mp.mpf(i) / 40 for i in range(-39, 1)]                 # (−1, 0]
kneser("sexp", "exp(1)", grund, dps=40)                          # Zwischenspeicher füllen


def T_kneser(x):
    k = int(mp.floor(x))
    r = x - k
    y = mp.re(kneser("sexp", "exp(1)", [r - 1], dps=40)[0]) if r != 0 else mp.mpf(0)
    return _schritte(y, k + 1)


def _schritte(y, m):
    """m-mal exp (m > 0) bzw. |m|-mal log (m < 0) — T(x+1) = e^{T(x)} in beide Richtungen."""
    for _ in range(max(m, 0)):
        y = mp.exp(y)
    for _ in range(max(-m, 0)):
        y = mp.log(y) if y > 0 else mp.mpf(-50)
    return y


def T_linear(x):
    k = int(mp.floor(x))
    y = x - k                                                      # auf [−1, 0): T = 1 + x, also T(r − 1) = r
    return _schritte(y, k + 1)


fig, ax = plt.subplots(1, 2, figsize=(9.6, 4))
ax[0].plot([float(x) for x in xs], [float(T_kneser(x)) for x in xs], color=BLAU, label="Kneser (glatt, analytisch)")
ax[0].plot([float(x) for x in xs], [float(T_linear(x)) for x in xs], "--", color=ROT, lw=1, label="„lineare“ Tetration (Knicke)")
ax[0].plot([-1, 0, 1, 2], [0, 1, E, E ** E], "o", color="#222", ms=4, label="ganze Höhen: 0, 1, e, e^e")
ax[0].set_ylim(-3, 22); ax[0].set_xlim(-1.95, 2.35); ax[0].axhline(0, color=GRAU, lw=0.5)
ax[0].set_xlabel("Höhe x"); ax[0].set_title("Basis e: zwei Antworten auf „halbe Höhen“"); ax[0].legend(fontsize=8, loc="upper left")
ein = ax[0].inset_axes([0.47, 0.36, 0.36, 0.3])
xd = [mp.mpf(i) / 40 for i in range(-40, 41)]
ein.plot([float(x) for x in xd], [float(T_linear(x) - T_kneser(x)) for x in xd], color=ROT, lw=1)
ein.axhline(0, color=GRAU, lw=0.5)
ein.set_title("Differenz linear − Kneser", fontsize=7); ein.tick_params(labelsize=6)
R = Regulaer(mp.sqrt(2), N=80, tol=mp.mpf("1e-3"))
xs2 = [mp.mpf(i) / 10 for i in range(-19, 81)]
ax[1].plot([float(x) for x in xs2], [float(R.tet(x)) for x in xs2], color=GRUEN, label="reguläre Tetration zur Basis √2")
ax[1].axhline(2, ls=":", color=GRAU)
ax[1].text(5.2, 2.05, "Fixpunkt 2", fontsize=8, color=GRAU)
ax[1].plot(range(0, 9), [float(R.tet(mp.mpf(k))) for k in range(0, 9)], "o", color="#222", ms=3, label="ganze Höhen")
ax[1].set_ylim(-3, 2.4); ax[1].set_xlabel("Höhe x"); ax[1].set_title("Basis √2: die Kurve kommt bei 2 zur Ruhe")
ax[1].legend(fontsize=8, loc="lower right")
speichern(fig, "k0_halbe_hoehe")

# ------------------------------------------------------------------ gr_young: Flächen unter T_e und neben slog_e
xs = [mp.mpf(i) / 60 for i in range(0, 61)]
ys = [float(T_kneser(x)) for x in xs]
fig, ax = plt.subplots(figsize=(5.2, 4.2))
xf = [float(x) for x in xs]
ax.fill_between(xf, 0, ys, color=BLAU, alpha=0.25, label="∫₀¹ T_e(x) dx ≈ 1.7145")
ax.fill_betweenx(ys, 0, xf, color=ROT, alpha=0.25, label="∫₁ᵉ slog_e(y) dy ≈ 1.0038")
ax.plot(xf, ys, color="#222", lw=1.4, label="y = T_e(x)")
ax.plot([0, 1, 1, 0, 0], [0, 0, E, E, 0], color=GRAU, lw=0.8, ls="--")
ax.text(0.33, 2.85, "Rechteck 1 × e", fontsize=9, color=GRAU)
ax.set_xlim(-0.05, 1.15); ax.set_ylim(0, 3.0); ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_title("Zwei Flächen ergeben zusammen genau e"); ax.legend(fontsize=8, loc="upper left")
speichern(fig, "gr_young")

# ================================================================== Etappe 2
from fixpunkt import AmFixpunkt  # noqa: E402

# ------------------------------------------------------------------ k4_koenigs: am Fixpunkt ist jeder Schritt eine Streckung
F2 = AmFixpunkt(mp.log(2) / 2, mp.mpf(2), N=60)
s1 = mp.re(F2.sigma(mp.mpf(1)))
lam = mp.log(2)
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.8))
z, bahn = mp.mpf(1), []
for n in range(0, 26):
    bahn.append(z)
    z = mp.sqrt(2) ** z
ax[0].semilogy(range(26), [float(abs(b - 2)) for b in bahn], "o", ms=4, color=ROT, label="Abstand |zₙ − 2| des Turms")
ax[0].semilogy(range(26), [float(abs(s1) * lam ** n) for n in range(26)], "-", color=BLAU, lw=1, label="|σ(1)|·λⁿ,  λ = ln 2")
ax[0].set_xlabel("Stockwerk n"); ax[0].set_title("Nahe am Fixpunkt: jeder Schritt ≈ mal λ"); ax[0].legend(fontsize=8)
xs = [mp.mpf(i) / 20 for i in range(0, 161)]
ax[1].plot([float(x) for x in xs], [float(abs(s1) * lam ** x) for x in xs], color=BLAU, lw=1.2, label="σ(T(x)) = σ(1)·λˣ (Betrag)")
ax[1].plot(range(9), [float(abs(F2.sigma(bahn[n]))) for n in range(9)], "o", color=ROT, ms=5, label="ganze Höhen (Turm)")
ax[1].plot([k + 0.5 for k in range(8)], [float(abs(s1) * lam ** (k + mp.mpf("0.5"))) for k in range(8)], "s", color=GRUEN, ms=4,
           label="halbe Höhen: einfach λ^{1/2} dazwischen")
ax[1].set_yscale("log"); ax[1].set_xlabel("Höhe x"); ax[1].set_title("In der Koenigs-Koordinate σ ist ein halber Schritt leicht")
ax[1].legend(fontsize=8)
speichern(fig, "k4_koenigs")

# ------------------------------------------------------------------ k4_blume: parabolischer Fixpunkt von e^u − 1
fig, ax = plt.subplots(figsize=(6.2, 4.4))
for w in np.linspace(0, 2 * np.pi, 48, endpoint=False):
    u = 0.45 * np.exp(1j * w)
    weg = [u]
    for _ in range(400):
        u = np.exp(u) - 1
        if abs(u) > 3:
            break
        weg.append(u)
    weg = np.array(weg)
    if abs(weg[-1]) < 0.05:
        ax.plot(weg.real, weg.imag, lw=0.7, color=BLAU, alpha=0.8)
        ax.plot([weg[0].real], [weg[0].imag], ".", color=BLAU, ms=3)
    else:
        ax.plot(weg[:6].real, weg[:6].imag, lw=0.7, color=ROT, alpha=0.7)
        ax.plot([weg[0].real], [weg[0].imag], ".", color=ROT, ms=3)
ax.plot([0], [0], "o", color="#222", ms=5)
ax.annotate("Fixpunkt 0, Multiplikator 1", (0, 0), xytext=(0.1, 0.32), fontsize=8, arrowprops=dict(arrowstyle="-", color=GRAU))
ax.set_aspect("equal"); ax.set_xlim(-0.6, 0.75); ax.set_ylim(-0.5, 0.5)
ax.set_title("e^u − 1: Bahnen kriechen von links an 0 heran (Blütenblatt)")
ax.plot([], [], color=BLAU, label="läuft in den Fixpunkt (anziehendes Blatt)")
ax.plot([], [], color=ROT, label="läuft weg (abstossendes Blatt)")
ax.legend(fontsize=8, loc="lower right")
speichern(fig, "k4_blume")

# ------------------------------------------------------------------ k5_vier: vier reguläre Superexponentiale zur Basis √2
F4 = AmFixpunkt(mp.log(2) / 2, mp.mpf(4), N=60)
fig, ax = plt.subplots(figsize=(7.4, 4.4))
xa = [mp.mpf(i) / 10 for i in range(-19, 81)]
ax.plot([float(x) for x in xa], [float(R.tet(x)) for x in xa], color=GRUEN, lw=1.6, label="am Fixpunkt 2, von unten: die Tetration (T(0) = 1)")
xb = [mp.mpf(i) / 10 for i in range(-60, 81)]
ax.plot([float(x) for x in xb], [float(mp.re(F2.Psi(lam ** x))) for x in xb], color=BLAU, lw=1.6, label="am Fixpunkt 2, von oben (zwischen 2 und 4)")
y0 = mp.re(F2.Psi(mp.mpf(1)))                                     # gleiche Normierung: beide Kurven bei x = 0 auf derselben Höhe
x0 = mp.findroot(lambda x: mp.re(F4.Psi(-(F4.lam) ** x)) - y0, mp.mpf(0))
ax.plot([float(x) for x in xb], [float(mp.re(F4.Psi(-(F4.lam) ** (x + x0)))) for x in xb], "--", color=ROT, lw=1.2,
        label="am Fixpunkt 4 (zwischen 2 und 4) — dieselbe Kurve bis auf 4·10⁻²⁵")
xc = [mp.mpf(i) / 20 for i in range(-120, 61)]
yc = [mp.re(F4.Psi((F4.lam) ** x)) for x in xc]
ax.plot([float(x) for x, y in zip(xc, yc) if y < 9], [float(y) for y in yc if y < 9], color=VIOLETT, lw=1.6, label="am Fixpunkt 4, nach oben (ganze Funktion)")
for y, t in ((2, "2"), (4, "4")):
    ax.axhline(y, ls=":", color=GRAU, lw=0.8)
ax.set_ylim(-2.5, 9); ax.set_xlim(-6, 8); ax.set_xlabel("x"); ax.set_title("Basis √2: vier Lösungen von F(x+1) = √2^F(x)")
ax.legend(fontsize=7.5, loc="upper left")
speichern(fig, "k5_vier")

# ================================================================== Etappe 3
from gemeinsam import gp_json  # noqa: E402

mp.mp.dps = 40
NK = 44


def _cheb(werte):
    k = [2 * mp.fsum(werte[j] * mp.cos(mp.pi * m * (j + mp.mpf(1) / 2) / NK) for j in range(NK)) / NK for m in range(NK)]
    k[0] /= 2
    return k


def _cheb_wert(k, t):
    b1 = b2 = mp.mpf(0)
    for c in reversed(k[1:]):
        b1, b2 = 2 * t * b1 - b2 + c, b1
    return t * b1 - b2 + k[0]


_kn = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / NK) / 2 for j in range(NK)]
_kA = _cheb([mp.re(w) for w in kneser("slog", "exp(1)", [mp.mpf(1) / 2 + t for t in _kn], dps=40)])   # slog_e auf [0, 1]
_kT2 = _cheb([mp.re(w) for w in kneser("sexp", "2", [-mp.mpf(1) / 2 + t for t in _kn], dps=40)])    # T_2 auf [−1, 0]


def slog_e(y):
    y, n = mp.mpf(y), 0
    while y > 1:
        y = mp.log(y); n += 1
    while y < 0:
        y = mp.exp(y); n -= 1
    return n + _cheb_wert(_kA, 2 * y - 1)


def T_2(x):
    k = int(mp.floor(x))
    y = _cheb_wert(_kT2, 2 * (x - k) - 1)
    for _ in range(k + 1):
        y = mp.power(2, y)
    return y


def Phi_e2(n, th):
    """slog_e(T_2(n+θ)) − (n+θ); ab n = 5 über die Logarithmenleiter (ln^m T_2(4+m+θ) = T_2(4+θ)·ln2 + ln ln2 für m ≥ 2)."""
    if n <= 4:
        return slog_e(T_2(n + th)) - n - th
    u = T_2(4 + th) * mp.log(2) + (mp.log(mp.log(2)) if n >= 6 else 0)
    return slog_e(u) - 4 - th


# ------------------------------------------------------------------ k6_regulaer_an_L: reguläre Superfunktion am komplexen Fixpunkt L_e
Le = mp.conj(-mp.lambertw(-1))
Le = Le if mp.im(Le) > 0 else mp.conj(Le)
FL = AmFixpunkt(1, Le, N=80)
s2 = FL.sigma(mp.mpf(2))                                           # F(0) = 2, F(−1) = ln 2, F(−2) = ln ln 2
x0 = slog_e(2)                                                     # Kneser erreicht 2 in der Höhe x0


def F_reg(x):
    return FL.Psi(FL.lam ** x * s2)


fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.2))
xr = [mp.mpf(i) / 40 for i in range(-400, 31)]
zr = [F_reg(x) for x in xr]
ax[0].plot([float(mp.re(z)) for z in zr], [float(mp.im(z)) for z in zr], color=ROT, lw=1.2, label="reguläre Lösung an L (läuft rückwärts in L)")
ax[0].plot([float(mp.re(z)) for z in zr], [-float(mp.im(z)) for z in zr], "--", color=ROT, lw=0.9, alpha=0.6, label="gespiegelt: dieselbe an L̄")
ax[0].plot([-3, 3.4], [0, 0], color=BLAU, lw=3, alpha=0.6, label="Kneser: bleibt auf der reellen Achse")
for x, t in ((0, "2"), (-1, "ln 2"), (-2, "ln ln 2")):
    z = F_reg(mp.mpf(x))
    ax[0].plot([float(mp.re(z))], [0], "o", color="#222", ms=4)
    ax[0].annotate(t, (float(mp.re(z)), 0), xytext=(float(mp.re(z)) - 0.1, -0.35), fontsize=8)
for z, t in ((Le, "L"), (mp.conj(Le), "L̄")):
    ax[0].plot([float(mp.re(z))], [float(mp.im(z))], "*", color="#222", ms=9)
    ax[0].annotate(t, (float(mp.re(z)), float(mp.im(z))), xytext=(float(mp.re(z)) + 0.12, float(mp.im(z)) + 0.1), fontsize=9)
ax[0].set_xlim(-3, 3.4); ax[0].set_ylim(-3.6, 3.6); ax[0].set_aspect("equal")
ax[0].set_xlabel("Re"); ax[0].set_ylabel("Im"); ax[0].set_title("Basis e in der komplexen Ebene"); ax[0].legend(fontsize=7, loc="lower left")
hs = [mp.mpf(i) / 20 for i in range(-38, 31)]
ax[1].plot([float(h) for h in hs], [float(T_kneser(h)) for h in hs], color=BLAU, lw=1.6, label="Kneser T_e(h) (reell)")
fh = [F_reg(h - x0) for h in hs]
ax[1].plot([float(h) for h in hs], [float(mp.re(z)) for z in fh], color=ROT, lw=1.2, label="regulär an L: Realteil")
ax[1].plot([float(h) for h in hs], [float(mp.im(z)) for z in fh], ":", color=ROT, lw=1.4, label="regulär an L: Imaginärteil")
for k in range(0, 3):
    ax[1].plot([float(x0 - k)], [float(T_kneser(x0 - k))], "o", color="#222", ms=4)
ax[1].axhline(0, color=GRAU, lw=0.5)
ax[1].set_ylim(-3, 6); ax[1].set_xlabel("Höhe h"); ax[1].set_title("Gleiche Werte nur bei 2, ln 2, ln ln 2")
ax[1].legend(fontsize=7.5, loc="upper left")
speichern(fig, "k6_regulaer_an_L")

# ------------------------------------------------------------------ k7_phase: slog_e(T_2(n+θ)) − (n+θ) läuft gegen μ + P(θ)
pm = gp_json("phi_modes.json")["bases"]


def grenzphase(b, th):
    s = mp.mpf(pm[b]["mu"])
    for k, a in enumerate(pm[b]["modes"], start=1):
        s += mp.re(mp.mpc(mp.mpf(a["re"]), mp.mpf(a["im"])) * mp.expjpi(2 * k * th))
    return s


fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.0))
ths = [mp.mpf(i) / 40 for i in range(41)]
farben = ["#c9ced8", "#a9b3c4", "#8291a8", "#56698a", "#2e4469"]
for n in range(5):
    ax[0].plot([float(t) for t in ths], [float(Phi_e2(n, t)) for t in ths], color=farben[n], lw=1.2, label=f"n = {n}")
ax[0].plot([float(t) for t in ths], [float(grenzphase("2", t)) for t in ths], color=ROT, lw=1.8, label="Grenzwert μ + P(θ)")
ax[0].set_xlabel("Nachkommateil θ der Höhe"); ax[0].set_ylabel("slog_e(T₂(n+θ)) − (n+θ)")
ax[0].set_ylim(-1.25, 0.45)
ax[0].set_title("Basis 2 in Höhen der Basis e übersetzt"); ax[0].legend(fontsize=7, loc="upper right", ncol=3)
th2 = [mp.mpf(i) / 80 for i in range(161)]
mu2, mu3 = mp.mpf(pm["2"]["mu"]), mp.mpf(pm["3"]["mu"])
ax[1].plot([float(t) for t in th2], [1e4 * float(grenzphase("2", t) - mu2) for t in th2], color=ROT, lw=1.5, label="P_e,2 (μ = −1.12840…)")
ax[1].plot([float(t) for t in th2], [1e4 * float(grenzphase("3", t) - mu3) for t in th2], color=GRUEN, lw=1.5, label="P_e,3 (μ = +0.19242…)")
tp = [mp.mpf(i) / 8 for i in range(17)]
ax[1].plot([float(t) for t in tp], [1e4 * float(Phi_e2(6, t - int(t)) - mu2) for t in tp], "o", color="#222", ms=3.5,
           label="nachgerechnet (n = 6), Abstand < 10⁻²³")
ax[1].axhline(0, color=GRAU, lw=0.5)
ax[1].set_xlabel("θ"); ax[1].set_ylabel("Welle P(θ) in Einheiten von 10⁻⁴"); ax[1].set_title("Was bleibt: eine winzige 1-periodische Welle")
ax[1].legend(fontsize=7.5, loc="lower left")
speichern(fig, "k7_phase")

# ------------------------------------------------------------------ k7_hub: die Hub-Spalte m(b) = μ_(e,b)
hub = gp_json("mu_hub.json")["values"]
bs = sorted(hub, key=float)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot([float(b) for b in bs], [float(hub[b]) for b in bs], "o-", color=BLAU, ms=3.5, lw=1.2, label="m(b) = μ_(e,b) aus der Tabelle")
ax.axvline(ETA, color=GRAU, ls=":", lw=0.8)
ax.text(ETA * 1.02, -11.3, "η: hier fällt m(b) nach −∞", fontsize=8, color=GRAU)
ax.axhline(0, color=GRAU, lw=0.5)
for b, t, dx, dy in (("2", "m(2) = −1.12840", 0.5, -2.6), ("3", "m(3) = +0.19242", 1.5, -1.6)):
    ax.plot([float(b)], [float(hub[b])], "o", color=ROT, ms=5)
    ax.annotate(t, (float(b), float(hub[b])), xytext=(float(b) + dx, float(hub[b]) + dy), fontsize=8, arrowprops=dict(arrowstyle="-", color=GRAU))
ax.plot([E], [0], "s", color="#222", ms=4)
ax.annotate("b = e: m = 0", (E, 0), xytext=(E * 2.2, 0.9), fontsize=8, arrowprops=dict(arrowstyle="-", color=GRAU))
ax.set_xscale("log"); ax.set_ylim(-12, 2.6); ax.set_xlabel("Basis b (logarithmisch)"); ax.set_ylabel("mittlere Höhenverschiebung gegen Basis e")
ax.set_title("Eine Spalte trägt alle Mittel: μ_(b,d) ≈ m(d) − m(b)"); ax.legend(fontsize=8, loc="lower right")
speichern(fig, "k7_hub")

# ================================================================== Etappe 4
# ------------------------------------------------------------------ k8_periode2: die Blase unter e^{−e} und die Rechteckwelle
re_ = np.linspace(-0.1, 0.085, 520)
im_ = np.linspace(-0.085, 0.085, 480)
B = re_[None, :] + 1j * im_[:, None]
C = np.log(B.astype(complex))
Z = np.ones_like(B)
with np.errstate(all="ignore"):
    for _ in range(3000):
        Z = np.exp(C * Z)
    Z1 = np.exp(C * Z)
    Z2 = np.exp(C * Z1)
    p2 = (np.abs(Z2 - Z) < 1e-9) & (np.abs(Z1 - Z) > 1e-7)
fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.3))
ax[0].contourf(re_, im_, p2.astype(float), levels=[0.5, 1.5], colors=[BLAU], alpha=0.55)
ax[0].contour(re_, im_, p2.astype(float), levels=[0.5], colors=[BLAU], linewidths=1)
w = np.linspace(0, 2 * np.pi, 400)
ax[0].plot(-0.0067418921 + 0.0722455899 * np.cos(w), 0.0722455899 * np.sin(w), "--", color=ROT, lw=1, label="bester Kreis (bis 2.13 % daneben)")
EE = math.exp(-E)
ax[0].plot([EE], [0], "s", color="#222", ms=4)
ax[0].annotate("e^(−e) = 0.0660", (EE, 0), xytext=(0.03, 0.018), fontsize=8, arrowprops=dict(arrowstyle="-", color=GRAU))
ax[0].plot([-0.0805281576], [0], "o", color="#222", ms=4)
ax[0].annotate("−0.0805", (-0.0805281576, 0), xytext=(-0.098, 0.016), fontsize=8)
ax[0].plot([0], [0], "x", color="#222", ms=5)
ax[0].set_aspect("equal"); ax[0].set_xlabel("Re b"); ax[0].set_ylabel("Im b")
ax[0].set_title("Basen mit Zweierwechsel: fast ein Kreis"); ax[0].legend(fontsize=7.5, loc="lower left")
bq = mp.mpf("0.03")


def T_p2(x):
    k = int(mp.floor(x))
    y = 1 - (1 - bq) * (x - k)                                    # Start: Gerade von T(0) = 1 nach T(1) = b
    for _ in range(k):
        y = bq ** y
    return y


xs8 = [mp.mpf(i) / 200 for i in range(0, 401)]
for n, f in zip((0, 1, 3, 10), ("#c9ced8", "#8291a8", "#56698a", "#2e4469")):
    ax[1].plot([float(x) for x in xs8], [float(T_p2(x + 2 * n)) for x in xs8], color=f, lw=1.2, label=f"T(x + {2 * n})")
for y, t in ((0.8213273730, "z_A"), (0.0561329673, "z_B"), (0.3226191538, "L (abstossend)")):
    ax[1].axhline(y, ls=":", color=GRAU, lw=0.8)
    ax[1].text(2.02, y, t, fontsize=7.5, va="center", color=GRAU)
ax[1].axvline(0.698330769, color=ROT, lw=0.8, ls="--"); ax[1].axvline(1.698330769, color=ROT, lw=0.8, ls="--")
ax[1].set_xlim(0, 2.35); ax[1].set_xlabel("x"); ax[1].set_title("b = 0.03: gerade Höhen werden zum Rechteck")
ax[1].legend(fontsize=7.5, loc="center right")
speichern(fig, "k8_periode2")
