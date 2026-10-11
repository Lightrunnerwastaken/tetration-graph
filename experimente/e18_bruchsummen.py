"""exp.bruchsummen — Basis √2: Grenzwertformel, Bruchsummen (Müller–Schleicher), fraktionale Produktformel.

 (1) T(x) = lim_n f^{−n}( a + λ^x (f^n(1) − a) )        (lineare Interpolation am Fixpunkt, zurückgezogen)
 (2) S(x) := Σ_{k=0}^{x−1} T(k) = x·a + Σ_{k≥0} [(T(k) − a) − (T(k+x) − a)]  ⇒  S(x+1) − S(x) = T(x), S(1) = 1
 (3) T'(x+h)/T'(x) = (ln b)^h · Π_{k=1}^{h} T(x+k) mit Bruchprodukt Π = exp(h ln a + Σ_{k≥1}[ln T(x+k) − ln T(x+k+h)])
     — gilt für die reguläre Lösung, verletzt für T̃ = T∘(id+θ) (Gegenprobe).
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 50
lauf = Lauf("exp.bruchsummen", __file__)
R = Regulaer(mp.sqrt(2), N=80, tol=mp.mpf("1e-3"))
a, lam, lb = R.a, R.lam_a, R.lb
f, finv = R.f, R.h


def T(x):
    return R.tet(x)


def Tstrich(x):
    return 1 / mp.diff(R.slog, T(x))


# (1) Grenzwertformel
fehler = {}
for n in (20, 60, 120, 200):
    y = mp.mpf(1)
    for _ in range(n):
        y = f(y)
    z = a + lam ** mp.mpf("0.5") * (y - a)
    for _ in range(n):
        z = finv(z)
    fehler[n] = abs(z - T(mp.mpf("0.5")))
lauf.pruefe("Grenzwertformel: Fehler bei n = 120 unter 1e-20", fehler[120] < mp.mpf(10) ** -20, {k: mp.nstr(v, 3) for k, v in fehler.items()})
# Rundungsverstärkung: Zurückziehen mit log_b verstärkt um 1/λ je Schritt ⇒ Gesamtfehler ≈ λ^n + ε·λ^{−n}, Optimum n ≈ dps·ln10/(2|ln λ|)
lauf.pruefe("Rundungsverstärkung sichtbar: bei n = 200 grösser als bei n = 120 (Optimum dazwischen)", fehler[200] > fehler[120],
            f"n_opt ≈ {mp.nstr(mp.mp.dps * mp.log(10) / (2 * abs(mp.log(lam))), 4)}")
quot = (fehler[120] / fehler[60]) ** (mp.mpf(1) / 60)
lauf.pruefe("Fehlerrate pro Schritt ≈ λ = ln 2", abs(quot - lam) < 0.02, mp.nstr(quot, 6))

# (2) Bruchsumme
K = 220


def S(x):
    tx = T(x)
    werte_x, werte_0 = [], []
    y0, yx = mp.mpf(1), tx
    for _ in range(K):
        werte_0.append(y0 - a); werte_x.append(yx - a)
        y0, yx = f(y0), f(yx)
    return x * a + mp.fsum(p - q for p, q in zip(werte_0, werte_x))


s_h = S(mp.mpf("0.5")); s_h1 = S(mp.mpf("1.5"))
lauf.pruefe("Bruchsumme: S(1/2 + 1) − S(1/2) = T(1/2)", abs(s_h1 - s_h - T(mp.mpf("0.5"))) < mp.mpf(10) ** -30, mp.nstr(abs(s_h1 - s_h - T(mp.mpf("0.5"))), 3))
lauf.pruefe("Bruchsumme: S(1) = T(0) = 1, S(2) = 1 + √2", abs(S(mp.mpf(1)) - 1) < 1e-30 and abs(S(mp.mpf(2)) - 1 - mp.sqrt(2)) < 1e-30)


# (3) Bruchprodukt-Formel
def bruchlog(Tfun, x, h):
    """Σ_{k=1}^{h} ln T(x+k) nach Müller–Schleicher (Grenzwert ln a)."""
    y1, yh = Tfun(x + 1), Tfun(x + 1 + h)
    s = mp.mpf(0)
    for _ in range(K):
        s += mp.log(y1) - mp.log(yh)
        y1, yh = f(y1), f(yh)
    return h * mp.log(a) + s


res = 0
for x in (mp.mpf(0), mp.mpf("-0.5")):
    for h in (mp.mpf("0.5"), mp.mpf("0.3"), mp.mpf("1.7")):
        links = mp.log(Tstrich(x + h)) - mp.log(Tstrich(x))
        rechts = h * mp.log(lb) + bruchlog(T, x, h)
        res = max(res, abs(links - rechts))
lauf.pruefe("reguläre Tetration: T'(x+h)/T'(x) = (ln b)^h · Bruchprodukt (h = 0.3, 0.5, 1.7)", res < mp.mpf(10) ** -20, mp.nstr(res, 3))

c = mp.mpf("0.05")
th = lambda x: c * (1 - mp.cos(2 * mp.pi * x)) / 2
Tt = lambda x: T(x + th(x))
Ttstrich = lambda x: Tstrich(x + th(x)) * (1 + c * mp.pi * mp.sin(2 * mp.pi * x))
gegen = []
for h in (mp.mpf("0.25"), mp.mpf("0.5")):
    x = mp.mpf("0.1")
    links = mp.log(Ttstrich(x + h)) - mp.log(Ttstrich(x))
    rechts = h * mp.log(lb) + bruchlog(Tt, x, h)
    gegen.append(abs(links - rechts))
lauf.pruefe("Gegenprobe T∘(id+θ), θ = 0.05·(1−cos 2πx)/2: Formel verletzt (Abweichung > 1e-3)", min(gegen) > 1e-3, [mp.nstr(g, 4) for g in gegen])

lauf.daten = {"grenzwert_fehler": {k: mp.nstr(v, 4) for k, v in fehler.items()}, "S_halb": mp.nstr(s_h, 25),
              "bruchprodukt_rest": mp.nstr(res, 3), "gegenprobe": [mp.nstr(g, 6) for g in gegen]}

xs = [mp.mpf(i) / 20 for i in range(0, 81)]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot([float(x) for x in xs], [float(S(x)) for x in xs], label="S(x) = Σ_{k<x} T(k) (Bruchsumme)")
teil = [mp.mpf(0)]
for k in range(4):
    teil.append(teil[-1] + T(mp.mpf(k)))
ax[0].plot(range(5), [float(t) for t in teil], "o", color="#b5523b", label="Teilsummen an ganzen Stellen")
ax[0].set_xlabel("x"); ax[0].set_title("Bruchsumme der √2-Tetration"); ax[0].legend(fontsize=8)
hs = [mp.mpf(i) / 20 for i in range(1, 41)]
x0 = mp.mpf("0.1")
d_reg = [float(abs(mp.log(Tstrich(x0 + h)) - mp.log(Tstrich(x0)) - h * mp.log(lb) - bruchlog(T, x0, h)) + 1e-40) for h in hs]
d_gest = [float(abs(mp.log(Ttstrich(x0 + h)) - mp.log(Ttstrich(x0)) - h * mp.log(lb) - bruchlog(Tt, x0, h)) + 1e-40) for h in hs]
ax[1].semilogy([float(h) for h in hs], d_reg, label="regulär")
ax[1].semilogy([float(h) for h in hs], d_gest, label="T∘(id+θ), θ klein")
ax[1].set_xlabel("h"); ax[1].set_title("Abweichung von der Bruchprodukt-Formel"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Grenzwertformel: Fehler {mp.nstr(fehler[120], 2)} bei n = 120, Rate ≈ λ; danach wächst die Rundung (Faktor 1/λ je Rückschritt). Bruchsumme S(1/2) = {mp.nstr(s_h, 16)}, "
                        f"S(x+1) − S(x) = T(x) exakt. Bruchprodukt-Formel regulär bis {mp.nstr(res, 2)}, für T∘(id+θ) verletzt um {mp.nstr(min(gegen), 3)}…{mp.nstr(max(gegen), 3)}.")
lauf.schreibe()
