"""exp.turmintegrale — I_n = ∫_0^1 ⁿx dx (Turm aus n Kopien von x), gerade/ungerade Grenzwerte, Lückenformel.

 Sophomore's Dream: I_2 = Σ (−1)^{k+1} k^{−k}.
 Grenzwerte (dominierte Konvergenz): lim I_{2n} = ∫_0^{e^{−e}} L₊ + ∫_{e^{−e}}^1 h,  lim I_{2n+1} = ∫_0^{e^{−e}} L₋ + ∫_{e^{−e}}^1 h,
 L₊ > L₋ der 2-Zyklus, h(x) = W(−ln x)/(−ln x).
 2-Zyklus parametrisiert durch u ∈ (0, 1/e): v ln v = u ln u, v = e^{W_0(u ln u)}, x = u^{1/v}.
 Andere Form: ∫_{e^{−e}}^1 h dx = 1 − e^{−e−1} − ∫_{1/e}^1 y^{1/y} dy.
 Vorhersage vor der Rechnung: I_{2n} − lim ≈ C n^{−3/2}.
"""
import mpmath as mp
import numpy as np

from gemeinsam import Lauf, plt

mp.mp.dps = 30
lauf = Lauf("exp.turmintegrale", __file__)
ee = mp.e ** (-mp.e)


def turm(x, n):
    t = x
    for _ in range(n - 1):
        t = x ** t
    return t


# Sophomore
I2 = mp.quad(lambda x: x ** x, [0, ee, 1])
reihe = mp.nsum(lambda k: (-1) ** (k + 1) / k ** k, [1, mp.inf])
lauf.pruefe("Sophomore's Dream: ∫_0^1 x^x dx = Σ (−1)^{k+1}/k^k", abs(I2 - reihe) < mp.mpf(10) ** -25, mp.nstr(I2, 20))
I_mp = {n: mp.quad(lambda x: turm(x, n), [0, ee, 0.3, 1]) for n in range(1, 7)}

# Grenzfunktionen
h_int_direkt = mp.quad(lambda x: mp.lambertw(-mp.log(x)).real / (-mp.log(x)), [ee, 0.5, 1])
h_int_form = 1 - mp.e ** (-mp.e - 1) - mp.quad(lambda y: y ** (1 / y), [1 / mp.e, 1])
lauf.pruefe("∫_{e^{−e}}^1 h = 1 − e^{−e−1} − ∫_{1/e}^1 y^{1/y} dy", abs(h_int_direkt - h_int_form) < mp.mpf(10) ** -22, mp.nstr(h_int_form, 20))


def zyklus(u):
    c = u * mp.log(u)
    v = mp.e ** mp.lambertw(c, 0).real
    return v, u ** (1 / v)


def dx_du(u):
    v, x = zyklus(u)
    nenner = mp.log(v) + 1
    if abs(nenner) < mp.mpf(10) ** (-mp.mp.dps + 5):      # Wurzelsingularität bei u = 1/e (integrierbar)
        return mp.mpf(0)
    vs = (mp.log(u) + 1) / nenner
    return x * (1 / (u * v) - mp.log(u) * vs / v ** 2)


oben = mp.quad(lambda u: zyklus(u)[0] * dx_du(u), [0, 0.1, 1 / mp.e])
unten = mp.quad(lambda u: u * dx_du(u), [0, 0.1, 1 / mp.e])
gerade, ungerade = oben + h_int_form, unten + h_int_form
luecke = oben - unten
# Gegencheck: 2-Zyklus direkt per Nullstellensuche bei einigen x
res = 0
for x in (mp.mpf("0.01"), mp.mpf("0.04"), mp.mpf("0.06")):
    hstern = mp.lambertw(-mp.log(x)).real / (-mp.log(x))      # abstossender Fixpunkt zwischen den Zyklus-Ästen
    lo = mp.findroot(lambda y: x ** (x ** y) - y, (mp.mpf("1e-15"), hstern * (1 - mp.mpf("1e-8"))), solver="illinois")
    hi = x ** lo
    uu = mp.findroot(lambda u: zyklus(u)[1] - x, (mp.mpf("1e-12"), 1 / mp.e - mp.mpf("1e-12")), solver="illinois")
    res = max(res, abs(uu - lo), abs(zyklus(uu)[0] - hi))
lauf.pruefe("Parametrisierung des 2-Zyklus = direkte Nullstellensuche", res < mp.mpf(10) ** -18, mp.nstr(res, 3))
bild2 = {"gerade": mp.mpf("0.676714273747536"), "ungerade": mp.mpf("0.631309370777713"), "luecke": mp.mpf("0.045404902969823")}
lauf.pruefe("Lücke stimmt mit Bild 2 auf 15 Stellen", abs(luecke - bild2["luecke"]) < 1e-14, mp.nstr(luecke, 18))
versatz_g, versatz_u = gerade - bild2["gerade"], ungerade - bild2["ungerade"]
lauf.pruefe("Bild 2: beide Grenzwerte um denselben Betrag zu klein (Fehler im gemeinsamen Teil ∫h)",
            abs(versatz_g - versatz_u) < 1e-14 and abs(versatz_g) > 1e-12, f"Versatz {mp.nstr(versatz_g, 4)} bzw. {mp.nstr(versatz_u, 4)}")

# Konvergenz von I_n (numpy, Gauss–Legendre in Paneelen, verdichtet bei 0 und bei e^{−e})
x_gl, w_gl = np.polynomial.legendre.leggauss(40)
ee_f = float(ee)
kanten = np.unique(np.concatenate([[0.0], np.geomspace(1e-14, 0.02, 40), np.linspace(0.02, ee_f, 30),
                                   ee_f + np.geomspace(1e-7, 0.02, 30), ee_f - np.geomspace(1e-7, 0.02, 30), np.linspace(ee_f + 0.02, 1, 40)]))
kanten = kanten[(kanten >= 0) & (kanten <= 1)]
a_, b_ = kanten[:-1], kanten[1:]
X = ((b_ - a_)[:, None] * (x_gl[None, :] + 1) / 2 + a_[:, None]).ravel()
Wt = ((b_ - a_)[:, None] * w_gl[None, :] / 2).ravel()
t = X.copy()
I_np = {1: float(np.sum(Wt * t))}
for n in range(2, 4097):
    t = X ** t
    I_np[n] = float(np.sum(Wt * t))
lauf.pruefe("numpy-Quadratur = mpmath für I_1…I_6 (1e-12)", max(abs(I_np[n] - float(I_mp[n])) for n in range(1, 7)) < 1e-12,
            f"{max(abs(I_np[n] - float(I_mp[n])) for n in range(1, 7)):.1e}")
ns = np.array([64, 128, 256, 512, 1024, 2048, 4096])
e_ger = np.array([I_np[n] - float(gerade) for n in ns])
e_ung = np.array([I_np[n + 1] - float(ungerade) for n in ns[:-1]] + [I_np[4095] - float(ungerade)])
p = -np.polyfit(np.log(ns), np.log(np.abs(e_ger)), 1)[0]
# Versuch verworfen: Richardson-Extrapolation der I_n bis n = 4096 erreicht nur ~2e-9 und kann den Versatz 2.4e-11 nicht entscheiden.
lauf.pruefe("Vorhersage: I_{2n} − lim fällt wie n^{−3/2} (Exponent 1.5 ± 0.1)", abs(p - 1.5) < 0.1, f"gemessen {p:.3f}")
lauf.daten = {"I": {n: mp.nstr(I_mp[n], 18) for n in I_mp}, "gerade": mp.nstr(gerade, 22), "ungerade": mp.nstr(ungerade, 22),
              "luecke": mp.nstr(luecke, 22), "int_h": mp.nstr(h_int_form, 22), "exponent": float(p),
              "fehler_gerade": {int(n): float(e) for n, e in zip(ns, e_ger)}}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
xs = np.linspace(1e-6, 1, 800)
for n in range(1, 7):
    tt = xs.copy()
    for _ in range(n - 1):
        tt = xs ** tt
    ax[0].plot(xs, tt, lw=1, label=f"ⁿx, n = {n}")
us = [mp.mpf(10) ** (-mp.mpf(k) / 10) * (1 / mp.e) for k in range(0, 140)]
zz = [zyklus(u) for u in us]
ax[0].plot([float(x) for _, x in zz], [float(v) for v, _ in zz], color="k", lw=2, label="2-Zyklus oben/unten")
ax[0].plot([float(x) for _, x in zz], [float(u) for u in us], color="k", lw=2)
ax[0].axvline(ee_f, ls=":", color="#888")
ax[0].set_xlabel("x"); ax[0].set_title("Türme auf (0,1) spalten bei e^{−e} in zwei Äste"); ax[0].legend(fontsize=7, ncol=2)
nn = np.arange(1, 41)
ax[1].plot(nn, [I_np[n] for n in nn], "o-", ms=3, label="I_n")
ax[1].axhline(float(gerade), ls="--", color="#b5523b", label="gerader Grenzwert")
ax[1].axhline(float(ungerade), ls="--", color="#2f8f83", label="ungerader Grenzwert")
ax[1].set_xlabel("n"); ax[1].set_title("I_n = ∫₀¹ ⁿx dx"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Sophomore bestätigt; gerade {mp.nstr(gerade, 16)}, ungerade {mp.nstr(ungerade, 16)}, Lücke {mp.nstr(luecke, 16)} "
                        f"(= Bild 2 auf 15 Stellen). ∫h = {mp.nstr(h_int_form, 14)}. I_n konvergiert wie n^(−{p:.3f}) (Vorhersage 3/2).")
lauf.schreibe()
