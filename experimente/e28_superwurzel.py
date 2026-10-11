"""exp.superwurzel — Superwurzeln: ssrt(x) = e^{W(ln x)}, Lösungszahl, Reihe um 1, dritte Superwurzel, unendliche Superwurzel x^{1/x}.

 ssrt(x): y^y = x  ⇔  y ln y = ln x  ⇔  ln y = W(ln x)  ⇒  y = e^{W(ln x)} = ln x / W(ln x).
 Reell: zwei Lösungen für e^{−1/e} < x < 1 (Zweige W₀, W₋₁), eine für x ≥ 1, keine für x < e^{−1/e}.
 Dritte Superwurzel: y^{y^y} = Y  ⇔  y = ssrt(ssrt(Y^y))  ⇒  Iteration y ← e^{W(W(y ln Y))} (Wikipedia).
 Unendliche Superwurzel: Turm von y^{1/y} konvergiert gegen y genau für 1/e ≤ y ≤ e (sonst gegen den anderen Fixpunkt oder gar nicht).
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 50
lauf = Lauf("exp.superwurzel", __file__)
E = mp.e


def ssrt(x, zweig=0):
    return mp.exp(mp.lambertw(mp.log(x), zweig))


# ------------------------------------------------------------------ Formel
xs = [mp.mpf(2), mp.mpf(10), mp.mpf(27), mp.mpf(1000), mp.mpf("0.9"), mp.exp(-1 / E) + mp.mpf("1e-3")]
res = max(abs(mp.re(ssrt(x)) ** mp.re(ssrt(x)) - x) for x in xs)
alt = max(abs(ssrt(x) - mp.log(x) / mp.lambertw(mp.log(x))) for x in xs if x != 1)
lauf.pruefe("ssrt(x) = e^{W(ln x)} löst y^y = x (6 Werte)", res < mp.mpf(10) ** -45, mp.nstr(res, 3))
lauf.pruefe("gleichwertige Form ln x / W(ln x)", alt < mp.mpf(10) ** -45, mp.nstr(alt, 3))
lauf.pruefe("ganze Fälle: ssrt(4) = 2, ssrt(27) = 3", abs(ssrt(4) - 2) + abs(ssrt(27) - 3) < mp.mpf(10) ** -45)

# ------------------------------------------------------------------ Anzahl reeller Lösungen
xmin = mp.exp(-1 / E)


def zaehle(x, n=4000):
    """Vorzeichenwechsel von y^y − x auf (0, 6]."""
    ys = [mp.mpf(i) / n * 6 for i in range(1, n + 1)]
    w = [y ** y - x for y in ys]
    return sum(1 for i in range(len(w) - 1) if w[i] * w[i + 1] < 0)


faelle = {"0.68 (unter e^{-1/e} = 0.6922)": (mp.mpf("0.68"), 0), "0.8": (mp.mpf("0.8"), 2), "0.95": (mp.mpf("0.95"), 2), "2": (mp.mpf(2), 1), "50": (mp.mpf(50), 1)}
gezaehlt = {k: zaehle(x) for k, (x, _) in faelle.items()}
lauf.pruefe("Lösungszahl: keine unter e^{−1/e}, zwei in (e^{−1/e}, 1), eine über 1",
            all(gezaehlt[k] == soll for k, (_, soll) in faelle.items()), gezaehlt)
x = mp.mpf("0.8")
y0, y1 = mp.re(ssrt(x, 0)), mp.re(ssrt(x, -1))
lauf.pruefe("die beiden Lösungen für x = 0.8 sind e^{W₀(ln x)} und e^{W₋₁(ln x)}, beide in (0, 1)",
            abs(y0 ** y0 - x) + abs(y1 ** y1 - x) < mp.mpf(10) ** -40 and 0 < y1 < 1 / E < y0 < 1, f"{mp.nstr(y1, 10)}, {mp.nstr(y0, 10)}")

# ------------------------------------------------------------------ Reihe um 1
koef = mp.taylor(lambda x: mp.re(ssrt(x)), 1, 7)
soll = [1, 1, -1, mp.mpf(3) / 2, mp.mpf(-17) / 6, mp.mpf(37) / 6, mp.mpf(-1759) / 120, mp.mpf(13279) / 360]
lauf.pruefe("Reihe um 1: 1 + (x−1) − (x−1)² + 3/2 (x−1)³ − 17/6 (x−1)⁴ + 37/6 − 1759/120 + 13279/360 …",
            max(abs(koef[i] - soll[i]) for i in range(8)) < mp.mpf(10) ** -30, [mp.nstr(c, 8) for c in koef])

# ------------------------------------------------------------------ dritte Superwurzel
dritte = {}
for Y in (2, 10, 100, 10 ** 6):
    Y = mp.mpf(Y)
    y = mp.mpf(1)
    for it in range(200):
        neu = mp.re(mp.exp(mp.lambertw(mp.lambertw(y * mp.log(Y)))))
        if abs(neu - y) < mp.mpf(10) ** -45:
            break
        y = neu
    dritte[int(Y)] = (y, it, abs(y ** (y ** y) - Y))
lauf.pruefe("dritte Superwurzel: Iteration y ← e^{W(W(y ln Y))} konvergiert und löst y^{y^y} = Y",
            all(r < mp.mpf(10) ** -38 * Y for Y, (_, _, r) in dritte.items()), {Y: f"{mp.nstr(v, 12)} ({i} Schritte)" for Y, (v, i, _) in dritte.items()})

# ------------------------------------------------------------------ unendliche Superwurzel
def turm(b, n=4000):
    t = mp.mpf(1)
    for _ in range(n):
        t = b ** t
    return t


innen = {str(y): turm(mp.mpf(y) ** (1 / mp.mpf(y))) for y in ("0.5", "1.5", "2", "2.6")}
lauf.pruefe("Turm von y^{1/y} konvergiert gegen y für 1/e ≤ y ≤ e (y = 0.5, 1.5, 2, 2.6)",
            all(abs(innen[k] - mp.mpf(k)) < mp.mpf(10) ** -8 for k in innen), {k: mp.nstr(v, 10) for k, v in innen.items()})
b3 = mp.mpf(3) ** (mp.mpf(1) / 3)
t3 = turm(b3)
a3 = mp.re(-mp.lambertw(-mp.log(b3)) / mp.log(b3))
lauf.pruefe("ausserhalb (y = 3): Turm von 3^{1/3} geht gegen den unteren Fixpunkt 2.4780…, nicht gegen 3",
            abs(t3 - a3) < mp.mpf(10) ** -20 and abs(t3 - 3) > mp.mpf("0.5"), mp.nstr(t3, 12))
ableitung_min = min(mp.diff(lambda y: y ** (1 / y), y) for y in mp.linspace(1 / E + mp.mpf("1e-6"), E - mp.mpf("1e-6"), 60))
lauf.pruefe("y ↦ y^{1/y} ist auf [1/e, e] streng wachsend (Umkehrung des unendlichen Turms)", ableitung_min > 0, mp.nstr(ableitung_min, 4))

lauf.daten = {"ssrt": {"2": mp.nstr(ssrt(2), 25), "3": mp.nstr(ssrt(3), 25), "10": mp.nstr(ssrt(10), 25)},
              "ssrt_0.8": [mp.nstr(y1, 20), mp.nstr(y0, 20)], "reihe_um_1": [mp.nstr(c, 15) for c in koef],
              "dritte": {Y: mp.nstr(v, 25) for Y, (v, _, _) in dritte.items()}, "turm_3_hoch_drittel": mp.nstr(t3, 20)}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
xg = [xmin + (mp.mpf(4) - xmin) * i / 400 for i in range(1, 401)]
ax[0].plot([float(x) for x in xg], [float(mp.re(ssrt(x))) for x in xg], color="#3d6fb6", label="ssrt(x) = e^{W₀(ln x)}")
xu = [x for x in xg if x < 1]
ax[0].plot([float(x) for x in xu], [float(mp.re(ssrt(x, -1))) for x in xu], "--", color="#b5523b", label="zweite Lösung e^{W₋₁(ln x)}")
ax[0].axvline(float(xmin), color="#8a93a3", lw=0.7, ls=":")
ax[0].text(float(xmin) + 0.03, 1.6, "e^{−1/e}", fontsize=8, color="#8a93a3")
ax[0].set_xlabel("x"); ax[0].set_ylabel("y mit y^y = x"); ax[0].set_title("Quadratische Superwurzel"); ax[0].legend(fontsize=8)
yg = [mp.mpf("0.2") + i * (mp.mpf("3.2") - mp.mpf("0.2")) / 300 for i in range(301)]
ax[1].plot([float(y) for y in yg], [float(y ** (1 / y)) for y in yg], color="#2f8f83", label="b = y^{1/y}")
ax[1].axvspan(float(1 / E), float(E), color="#2f8f83", alpha=0.08, label="1/e ≤ y ≤ e: Turm von b konvergiert gegen y")
ax[1].plot([float(E)], [float(E ** (1 / E))], "o", color="#222", ms=4)
ax[1].text(float(E) - 0.5, float(E ** (1 / E)) + 0.06, "(e, η)", fontsize=8)
ax[1].set_xlabel("y"); ax[1].set_ylabel("Basis b"); ax[1].set_title("Unendliche Superwurzel: b = y^{1/y}"); ax[1].legend(fontsize=7.5, loc="lower right")
lauf.bild(fig)
lauf.zusammenfassung = (f"ssrt(2) = {mp.nstr(ssrt(2), 12)}, ssrt(10) = {mp.nstr(ssrt(10), 12)}; Reihe um 1 bestätigt; Lösungszahl 0/2/1; "
                        f"dritte Superwurzel von 2: {mp.nstr(dritte[2][0], 12)}; unendliche Superwurzel = y^(1/y) genau auf [1/e, e].")
lauf.schreibe()
