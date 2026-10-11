"""exp.tetra_euler — η₄ und e₄: wo die Kneser-Tetration T_b die Gerade y = x berührt (eine Stufe über η und e).

 Für jede Basis b: Taylorreihe von T_b um 3.0885, Stelle x*(b) mit T_b'(x*) = 1 (Minimum von T_b(x) − x), g(b) = T_b(x*) − x*.
 η₄ ist die Nullstelle von g (Sekantenverfahren in b), e₄ = x*(η₄). Dieselbe Methode wie S. Levenstein (Tetration Forum tid 262,
 31.10.2010, 32 Stellen: 1.6353244967152763993453446183062 und 3.0885322718067176544821807826411).
 Drei Kennzeichnungen: (1) Berührung, (2) Maximum der Selbst-Tetrawurzel y ↦ b(y) mit T_b(y) = y, (3) b↑↑↑∞ konvergiert genau für b ≤ η₄.
 Fixpunktzahl: drei für η < b < η₄, einer für b > η₄. Unter η (regulär, Basis √2): Fixpunkt > 0 ist anziehend, √2↑↑↑∞ konvergiert.
"""
import mpmath as mp

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, kneser, plt
from stufen import basis_text, fixpunkt, p_abl, p_nullstelle, p_verschieben, p_wert, taylor

mp.mp.dps = 50
lauf = Lauf("exp.tetra_euler", __file__)
LEV = {"eta4": mp.mpf("1.6353244967152763993453446183062"), "e4": mp.mpf("3.0885322718067176544821807826411"),
       "L": mp.mpf("-1.6408725757165933485612321510790"), "lamL": mp.mpf("4.8060057543017516963843938970331")}
X0 = mp.mpf("3.0885")


def T(b, xs):
    return [mp.re(v) for v in kneser("sexp", basis_text(b), xs)]


def minimum(b):
    """x* mit T_b'(x*) = 1 nahe X0 und g = T_b(x*) − x*, dazu die Reihe um x*."""
    a = taylor(b, X0, mp.mpf("0.6"), N=96, K=45)
    da = p_abl(a)
    g1 = list(da)
    g1[0] -= 1
    t = p_nullstelle(g1)
    return X0 + t, p_wert(a, t) - (X0 + t), p_verschieben(a, t)


# ------------------------------------------------------------------ η₄ durch Sekanten in b
b0, b1 = mp.mpf("1.6353"), mp.mpf("1.6354")
g0, g1 = minimum(b0)[1], minimum(b1)[1]
for _ in range(12):
    b2 = b1 - g1 * (b1 - b0) / (g1 - g0)
    b0, g0 = b1, g1
    b1 = b2
    g1 = minimum(b1)[1]
    if abs(b1 - b0) < mp.mpf(10) ** -36:
        break
eta4 = b1
e4, g_e4, a4 = minimum(eta4)
print(f"  η₄ = {mp.nstr(eta4, 35)}\n  e₄ = {mp.nstr(e4, 35)}  (g = {mp.nstr(g_e4, 3)})")
lauf.pruefe("η₄ wie Levenstein (2010, 32 Stellen) auf 1e-28", abs(eta4 - LEV["eta4"]) < mp.mpf(10) ** -28, f"{mp.nstr(eta4, 32)} (Abstand {mp.nstr(eta4 - LEV['eta4'], 2)})")
lauf.pruefe("e₄ wie Levenstein auf 1e-26", abs(e4 - LEV["e4"]) < mp.mpf(10) ** -26, f"{mp.nstr(e4, 32)} (Abstand {mp.nstr(e4 - LEV['e4'], 2)})")
lauf.pruefe("Berührung: T(e₄) = e₄, T'(e₄) = 1, T''(e₄) > 0", abs(a4[0] - e4) < mp.mpf(10) ** -30 and abs(a4[1] - 1) < mp.mpf(10) ** -30 and a4[2] > 0,
            f"T''(e₄)/2 = {mp.nstr(a4[2], 12)}")

# negativer Fixpunkt bei η₄
Ln, an = fixpunkt(eta4, mp.mpf("-1.6409"), mp.mpf("0.2"), N=96, K=40)
lauf.pruefe("negativer Fixpunkt bei η₄: L = −1.64087257571659… mit Multiplikator 4.80600575430175… (Levenstein)",
            abs(Ln - LEV["L"]) < mp.mpf(10) ** -26 and abs(an[1] - LEV["lamL"]) < mp.mpf(10) ** -24, f"{mp.nstr(Ln, 25)}, λ = {mp.nstr(an[1], 20)}")

# ------------------------------------------------------------------ Fixpunktzahl: Vorzeichenwechsel von T_b(x) − x auf (−2, x_max]
def fixpunkte(b, schritt=mp.mpf("0.25")):
    """Zählt Vorzeichenwechsel von T_b(x) − x, bis T_b(x) > 60 (danach wächst T_b schneller als x, kein Schnitt mehr)."""
    x, vor, wechsel, stellen = mp.mpf("-1.99"), None, 0, []
    while True:
        xs = [x + k * schritt for k in range(8)]
        for xi, ti in zip(xs, T(b, xs)):
            d = ti - xi
            if vor is not None and d * vor < 0:
                wechsel += 1
                stellen.append(mp.nstr(xi - schritt / 2, 3))
            vor = d
            if ti > 60 and d > 0:
                return wechsel, stellen
        x = xs[-1] + schritt


zaehl = {}
for b in ("1.5", "1.6", "1.7", "2", "exp(1)"):
    zaehl[b] = fixpunkte(b)
lauf.pruefe("drei reelle Fixpunkte für η < b < η₄ (1.5, 1.6), einer für b > η₄ (1.7, 2, e)",
            [zaehl[b][0] for b in zaehl] == [3, 3, 1, 1, 1], {b: f"{z[0]} bei {z[1]}" for b, z in zaehl.items()})
g_unter, g_ueber = minimum(eta4 - mp.mpf("1e-4"))[1], minimum(eta4 + mp.mpf("1e-4"))[1]
lauf.pruefe("nahe η₄: Minimum von T_b(x) − x negativ für η₄ − 10⁻⁴ (zwei Schnitte bei e₄), positiv für η₄ + 10⁻⁴ (keiner)",
            g_unter < 0 < g_ueber, f"{mp.nstr(g_unter, 4)} / {mp.nstr(g_ueber, 4)}")

# ------------------------------------------------------------------ (2) Maximum der Selbst-Tetrawurzel
def selbstwurzel(y):
    """b mit T_b(y) = y (Sekanten in b)."""
    c0, c1 = mp.mpf("1.62"), mp.mpf("1.63")
    h0, h1 = T(c0, [y])[0] - y, T(c1, [y])[0] - y
    for _ in range(30):
        c2 = c1 - h1 * (c1 - c0) / (h1 - h0)
        c0, h0, c1 = c1, h1, c2
        h1 = T(c1, [y])[0] - y
        if abs(c1 - c0) < mp.mpf(10) ** -25:
            break
    return c1


sw = {mp.nstr(y, 4): selbstwurzel(y) for y in (mp.mpf("2.6"), mp.mpf("2.9"), mp.mpf("3.3"), mp.mpf("3.6"))}
lauf.pruefe("Selbst-Tetrawurzel b(y) (T_b(y) = y) ist kleiner als η₄ für y ≠ e₄ — Maximum η₄ bei e₄",
            all(v < eta4 for v in sw.values()), {k: mp.nstr(v, 12) for k, v in sw.items()})

# ------------------------------------------------------------------ (3) unendliche Pentation b↑↑↑n = T_b^n(1)
def turm5(b, n, grenze=5):
    """1, T(1), T(T(1)), … bis n Schritte oder über `grenze` (T_b(7.5) läuft im Motor über)."""
    folge = [mp.mpf(1)]
    for _ in range(n):
        folge.append(T(b, [folge[-1]])[0])
        if folge[-1] > grenze:
            break
    return folge


f16 = turm5("1.6", 40)
x1 = fixpunkt("1.6", f16[-1], mp.mpf("0.3"), N=64, K=30)
lauf.pruefe("b = 1.6 < η₄: 1.6↑↑↑n wächst monoton gegen den unteren positiven Fixpunkt (anziehend, T' < 1)",
            all(f16[i] < f16[i + 1] for i in range(len(f16) - 1)) and abs(f16[-1] - x1[0]) < mp.mpf(10) ** -10 and x1[1][1] < 1,
            f"Grenzwert {mp.nstr(x1[0], 20)}, T' = {mp.nstr(x1[1][1], 8)}, nach {len(f16) - 1} Schritten")
f17 = turm5("1.7", 60)
lauf.pruefe("b = 1.7 > η₄: 1.7↑↑↑n wächst an e₄ vorbei über 5 (danach T_b(x) > x ohne Schnitt → unbeschränkt)", f17[-1] > 5,
            f"{len(f17) - 1} Schritte: {[mp.nstr(v, 5) for v in f17[-4:]]}")

# unter η (reguläre Tetration, Basis √2): positiver Fixpunkt anziehend
u = mp.log(mp.sqrt(2))
F2 = AmFixpunkt(u, mp.mpf(2), N=80)
s1 = F2.sigma(mp.mpf(1))
Treg = lambda z: mp.re(F2.Psi(F2.lam ** z * s1))
xr = mp.findroot(lambda x: Treg(x) - x, mp.mpf("1.8"))
dr = mp.diff(Treg, xr)
f2 = [mp.mpf(1)]
for _ in range(60):
    f2.append(Treg(f2[-1]))
lauf.pruefe("unter η (Basis √2, regulär): positiver Fixpunkt von T anziehend, √2↑↑↑∞ konvergiert dorthin",
            dr < 1 and abs(f2[-1] - xr) < mp.mpf(10) ** -12, f"{mp.nstr(xr, 20)}, T' = {mp.nstr(dr, 10)}")

lauf.daten = {"eta4": mp.nstr(eta4, 35), "e4": mp.nstr(e4, 35), "T2_halbe": mp.nstr(a4[2], 20), "L_neg": mp.nstr(Ln, 30),
              "lam_neg": mp.nstr(an[1], 25), "fixpunkte": {b: z[0] for b, z in zaehl.items()},
              "selbstwurzel": {k: mp.nstr(v, 20) for k, v in sw.items()}, "pent_1.6": mp.nstr(x1[0], 25), "Tstrich_1.6": mp.nstr(x1[1][1], 15),
              "sqrt2_fix": mp.nstr(xr, 25), "sqrt2_Tstrich": mp.nstr(dr, 15)}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
gx = [mp.mpf("-1.95") + mp.mpf(i) / 20 for i in range(0, 125)]
for b, farbe, name in (("1.6", "#2f8f83", "b = 1.6"), (eta4, "#b5523b", "b = η₄"), ("1.7", "#7a49c2", "b = 1.7")):
    ys = T(b, gx)
    ax[0].plot([float(x) for x in gx], [float(y) for y in ys], color=farbe, label=name)
ax[0].plot([-2, 4.3], [-2, 4.3], ":", color="#8a93a3", label="y = x")
ax[0].plot([float(e4)], [float(e4)], "o", color="#b5523b", ms=4)
ax[0].set_ylim(-2.5, 6); ax[0].set_xlim(-2, 4.3); ax[0].set_xlabel("x"); ax[0].set_title("T_b berührt y = x bei e₄"); ax[0].legend(fontsize=8)
ys3 = [mp.mpf("2.2") + mp.mpf(i) / 50 for i in range(0, 86)]
ETA = mp.e ** (1 / mp.e)
ax[1].plot([float(y) for y in ys3], [float(y ** (1 / y) - ETA) for y in ys3], color="#3d6fb6", label="Stufe 3: y^{1/y} − η (Max. bei e)")
ax[1].plot([float(mp.e)], [0], "o", color="#3d6fb6", ms=4)
sw_bild = dict(sw)
for y in ("2.3", "2.45", "2.75", "3.45", "3.8"):
    sw_bild[y] = selbstwurzel(mp.mpf(y))
sx = sorted([(mp.mpf(k), v) for k, v in sw_bild.items()] + [(e4, eta4)])
ax[1].plot([float(p[0]) for p in sx], [float(p[1] - eta4) for p in sx], "o-", color="#b5523b", ms=3.5, label="Stufe 4: b(y) − η₄, T_b(y) = y (Max. bei e₄)")
ax[1].axhline(0, color="#8a93a3", lw=0.6)
ax[1].set_xlabel("y"); ax[1].set_ylabel("Abstand zur Grenzbasis"); ax[1].set_title("Eine Stufe höher: (e, η) → (e₄, η₄)"); ax[1].legend(fontsize=7.5, loc="lower center")
lauf.bild(fig)
lauf.zusammenfassung = (f"η₄ = {mp.nstr(eta4, 25)}, e₄ = {mp.nstr(e4, 25)} (Levenstein 2010 bestätigt); Fixpunktzahl 3 → 1 bei η₄; "
                        f"1.6↑↑↑∞ = {mp.nstr(x1[0], 12)}, 1.7↑↑↑n divergiert; Selbst-Tetrawurzel maximal bei e₄.")
lauf.schreibe()
