"""exp.turmintegral_periodisch — das Turm-Integral bei nicht ganzer Höhe: I(s) = ∫_0^1 T_x(s) dx.

 Für jede Basis x ∈ (0,1) braucht man eine stetige reelle Tetration T_x (T_x(0) = 1, T_x(s+1) = x^{T_x(s)}). Drei Familien:
   A „linear“:       T_x(s) = x^s auf [0,1]            (aus der linearen Tetration auf (−1,0])
   B „geradlinig“:   T_x(s) = 1 − s(1 − x) auf [0,1]
   C „C¹-glatt“:     T_x(s) = x^{φ(s)}, φ(s) = s + c s(1−s), c = (1 − ln x)/(1 + ln x)   (Steigung passt an den ganzen Stellen)
 Nach fp.satz.rechteckwelle geht T_x(2m + r) für x < e^{−e} gegen eine Rechteckwelle mit Sprung an der Kreuzungshöhe s*(x)
 (T_x(s*) = L(x), abstossender Fixpunkt), für x ≥ e^{−e} gegen h(x). Daraus: I(2m + r) → I_∞(r), 2-periodisch, mit
   I_∞(0) = lim I_2n, I_∞(1) = lim I_2n+1, I_∞(r+1) = 2M − I_∞(r), monoton auf [0,1], Form = Verteilung von s*(x).
"""
import mpmath as mp
import numpy as np

from gemeinsam import Lauf, plt

mp.mp.dps = 30
lauf = Lauf("exp.turmintegral_periodisch", __file__)
EE = float(mp.e ** (-mp.e))
I_GERADE = mp.mpf("0.67671427377181279519")
I_UNGERADE = mp.mpf("0.63130937080198819539")
M = (I_GERADE + I_UNGERADE) / 2

# 2-Zyklus über u ∈ (0, 1/e): v ln v = u ln u, x = u^(1/v); L₋ = u, L₊ = v; abstossender Fixpunkt L = W(−ln x)/(−ln x)
gl_u, gl_w = np.polynomial.legendre.leggauss(1500)


def zyklus(u):
    v = mp.exp(mp.lambertw(u * mp.log(u)).real)
    return v, mp.power(u, 1 / v)


# Knoten in u mit Substitution u = (1/e)·(1 − τ²) (dichter nahe 1/e, wo x → e^{−e})
tau = (gl_u + 1) / 2
us = [mp.mpf(1) / mp.e * (1 - mp.mpf(t) ** 2) for t in tau]
dudtau = [2 * mp.mpf(t) / mp.e for t in tau]
daten = []
for u, du, w in zip(us, dudtau, gl_w):
    v, x = zyklus(u)
    dxdu = mp.diff(lambda uu: zyklus(uu)[1], u)
    L = mp.lambertw(-mp.log(x)).real / (-mp.log(x))
    daten.append((float(x), float(v - u), float(-dxdu * du * w / 2), float(L)))   # u ↓ ⇒ x ↓: Vorzeichen so, dass Gewicht > 0
X = np.array([d[0] for d in daten]); GAP = np.array([d[1] for d in daten]); GEW = np.abs(np.array([d[2] for d in daten])); LF = np.array([d[3] for d in daten])
luecke = float((GAP * GEW).sum())
lauf.pruefe("Gewichte in u geben die bekannte Lücke ∫(L₊ − L₋)dx = 0.0454049…", abs(luecke - float(I_GERADE - I_UNGERADE)) < 1e-7, f"{luecke:.12f}")


def kreuzung(familie, x, L):
    lx = np.log(x)
    if familie == "A":
        return L.copy()                                  # x^s = L = x^L  ⇒  s* = L
    if familie == "B":
        return (1 - L) / (1 - x)
    c = (1 - lx) / (1 + lx)                              # φ(s) = L:  c s² − (1 + c) s + L = 0
    disk = np.sqrt((1 + c) ** 2 - 4 * c * L)
    s1, s2 = ((1 + c) - disk) / (2 * c), ((1 + c) + disk) / (2 * c)
    return np.where((s1 >= 0) & (s1 <= 1), s1, s2)


def I_grenz(familie, r):
    """I_∞(r) für r ∈ [0, 2]."""
    s = kreuzung(familie, X, LF)
    if r <= 1:
        return float(I_GERADE) - float((GAP * GEW)[s < r].sum())
    return 2 * float(M) - I_grenz(familie, r - 1)


rs = np.linspace(0, 2, 201)
kurven = {f: np.array([I_grenz(f, r) for r in rs]) for f in "ABC"}
ok_enden = all(abs(k[0] - float(I_GERADE)) < 1e-12 and abs(k[100] - float(I_UNGERADE)) < 1e-7 and abs(k[200] - float(I_GERADE)) < 1e-7 for k in kurven.values())
lauf.pruefe("bei jeder Interpolation: I_∞(0) = lim I₂ₙ, I_∞(1) = lim I₂ₙ₊₁, I_∞(2) = lim I₂ₙ", ok_enden,
            {f: [round(k[0], 9), round(k[100], 9)] for f, k in kurven.items()})
mono = all(np.all(np.diff(k[:101]) <= 1e-9) and np.all(np.diff(k[100:]) >= -1e-9) for k in kurven.values())   # 1e-9: Quadraturfehler der Lücke
lauf.pruefe("monoton: fallend auf [0,1], steigend auf [1,2]", mono)
# Fourier über eine Periode: Mittel M, gerade Harmonische 0
fourier = {}
for f, k in kurven.items():
    werte = k[:-1]
    c = np.fft.rfft(werte) / len(werte)
    fourier[f] = {"Mittel": round(float(c[0].real), 9), "|c2|": f"{abs(c[2]):.1e}", "|c1|": f"{abs(c[1]):.4e}", "|c3|": f"{abs(c[3]):.4e}"}
lauf.pruefe("Mittel = M = 0.654011822…, gerade Harmonische verschwinden (Halbperioden-Antisymmetrie)",
            all(abs(v["Mittel"] - float(M)) < 1e-6 and float(v["|c2|"]) < 1e-9 for v in fourier.values()), fourier)
mitte = {f: round(float(k[50]), 6) for f, k in kurven.items()}
lauf.pruefe("die Form hängt von der Interpolation ab: I_∞(½) verschieden für A, B, C", len(set(mitte.values())) == 3, mitte)
# A: Abfall nur auf [0, 1/e], dort geschlossen I_∞(r) = lim I₂ₙ − ∫_0^{r^{1/r}} (L₊ − L₋) dx
kA = kurven["A"]
plateau = np.all(np.abs(kA[(rs > 1 / np.e + 0.01) & (rs <= 1)] - float(I_UNGERADE)) < 1e-7)
r0 = 0.25
xr = r0 ** (1 / r0)
geschlossen = float(I_GERADE) - float((GAP * GEW)[X < xr].sum())
lauf.pruefe("lineare Interpolation: s*(x) = L(x) genau, Abfall nur für r < 1/e, I_∞(r) = lim I₂ₙ − ∫₀^{r^{1/r}}(L₊ − L₋)dx",
            plateau and abs(geschlossen - I_grenz("A", r0)) < 1e-12, {"I_∞(0.25)": round(geschlossen, 9)})


# Endliche Höhen: I(2m + r) mit Gauss–Legendre in x, gegen den Grenzwert
gx, gw = np.polynomial.legendre.leggauss(4000)


def knoten(a, b):
    return (b - a) / 2 * gx + (a + b) / 2, (b - a) / 2 * gw


x1, w1 = knoten(1e-12, EE)
x2, w2 = knoten(EE, 1.0)
XS, WS = np.concatenate([x1, x2]), np.concatenate([w1, w2])


def T_start(familie, x, r):
    """T_x(r) für r ∈ [0, 2)."""
    if r >= 1:
        return x ** T_start(familie, x, r - 1)
    if familie == "A":
        return x ** r
    if familie == "B":
        return 1 - r * (1 - x)
    lx = np.log(x)
    c = (1 - lx) / (1 + lx)
    return x ** (r + c * r * (1 - r))


def I_endlich(familie, s):
    m, r = divmod(s, 2)
    t = T_start(familie, XS, r)
    for _ in range(int(round(m)) * 2):
        t = XS ** t
    return float((t * WS).sum())


abst = {}
for m in (5, 20, 80):
    abst[m] = max(abs(I_endlich("A", 2 * m + r) - I_grenz("A", r)) for r in (0.1, 0.5, 1.3))
lauf.pruefe("I(2m + r) nähert sich I_∞(r) (lineare Familie, m = 5, 20, 80)", abst[80] < abst[5], {m: f"{v:.2e}" for m, v in abst.items()})

lauf.daten = {"M": mp.nstr(M, 20), "fourier": fourier, "I_halb": mitte}
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
namen = {"A": "linear: x^s", "B": "geradlinig: 1 − s(1−x)", "C": "C¹-glatt"}
for f, k in kurven.items():
    ax[0].plot(rs, k, label=namen[f])
ax[0].axhline(float(M), ls=":", color="#888")
ax[0].scatter([0, 1, 2], [float(I_GERADE), float(I_UNGERADE), float(I_GERADE)], color="k", zorder=5, s=18)
ax[0].set_xlabel("r (Höhe 2m + r, m → ∞)"); ax[0].set_title("Grenzkurve I_∞(r): gleiche Ecken, verschiedene Form")
ax[0].legend(fontsize=8)
for m in (2, 5, 20):
    ax[1].plot(rs, [I_endlich("A", 2 * m + r) for r in rs], label=f"m = {m}")
ax[1].plot(rs, kurven["A"], "k--", lw=1, label="Grenze")
ax[1].set_xlabel("r"); ax[1].set_title("lineare Familie: I(2m + r) für wachsendes m")
ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"I(2m + r) → 2-periodische Grenzkurve I_∞(r) mit festen Ecken lim I₂ₙ, lim I₂ₙ₊₁, Mittel M = {mp.nstr(M, 12)}, "
                        "nur ungeraden Harmonischen und monotonen Flanken; die Form ist die Verteilung der Kreuzungshöhen und hängt von der Interpolation ab.")
lauf.schreibe()
