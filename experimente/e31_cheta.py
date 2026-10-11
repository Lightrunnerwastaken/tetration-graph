"""exp.cheta — obere Superexponentialfunktionen: Basis √2 am oberen Fixpunkt 4, Basis η über e ("cheta", J. D. Fox 2009).

 √2: usexp(x) = Ψ₄(λ₄^x · c), λ₄ = 4 ln √2 = 2 ln 2 > 1 (abstossend), normiert usexp(0) = 5 (A + 1); ganz, Periode 2πi/ln λ₄.
 η: cheta(z) = exp_η^{∘z}(e²) mit cheta(−1) = 2e, cheta(0) = e², cheta(1) = e^e. Über die Kehrwert-Koordinate Z = 2e/(e − y) < 0 und
 die Fatou-Reihe Φ(Z) = Z − (1/3) ln|Z| + Σ c_k Z^{−k} (dieselben rationalen c_k wie unten, exp.eta_digamma):
 C_ch = lim_n [Φ(Z(cheta(−n))) + n], cheta(x) = e − 2e/Z mit Φ(Z) = x + C_ch (für x → −∞), sonst vorwärts iterieren.
"""
import mpmath as mp
import sympy as sp

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, plt

mp.mp.dps = 50
lauf = Lauf("exp.cheta", __file__)
E = mp.e
ETA = mp.e ** (1 / mp.e)

# ------------------------------------------------------------------ Basis √2, oberer Fixpunkt 4
b = mp.sqrt(2)
F4 = AmFixpunkt(mp.log(b), mp.mpf(4), N=80)
c4 = F4.sigma(mp.mpf(5))


def usexp(z):
    return F4.Psi(F4.lam ** z * c4)


werte = {"0": usexp(0), "1": usexp(1), "-1": usexp(-1), "1/2": usexp(mp.mpf(1) / 2)}
lauf.pruefe("√2 oben: usexp(0) = 5, usexp(1) = √2⁵, usexp(−1) = log_√2 5",
            abs(werte["0"] - 5) + abs(werte["1"] - b ** 5) + abs(werte["-1"] - mp.log(5) / mp.log(b)) < mp.mpf(10) ** -40,
            {k: mp.nstr(mp.re(v), 15) for k, v in werte.items()})
zk = mp.mpc("0.3", "2.0")
lauf.pruefe("√2 oben: Gleichung usexp(z+1) = √2^usexp(z) auch komplex (z = 0.3 + 2i)", abs(usexp(zk + 1) - b ** usexp(zk)) < mp.mpf(10) ** -38)
per = 2j * mp.pi / mp.log(mp.re(F4.lam))
lauf.pruefe("√2 oben: periodisch mit 2πi/ln λ₄ (Periode 19.24 i)", abs(usexp(zk + per) - usexp(zk)) < mp.mpf(10) ** -38, mp.nstr(mp.im(per), 6))
lauf.pruefe("√2 oben: usexp(x) → 4 für x → −∞ (Abstand bei x = −30 wie λ₄^x)", abs(usexp(-30) - 4) < 2 * abs(c4) * abs(F4.lam) ** -30,
            mp.nstr(usexp(-30) - 4, 4))

# ------------------------------------------------------------------ Basis η: Fatou-Reihe (rational)
K = 24
t = sp.symbols("t")
eps = sp.series(sp.coth(t) - 1 / t, t, 0, K + 3).removeO()
cs = sp.symbols(f"c1:{K + 1}")
G = eps - sp.Rational(1, 3) * sp.log(1 + t + t * eps)
for k in range(1, K + 1):
    G += cs[k - 1] * t ** k * ((1 + t + t * eps) ** (-k) - 1)
G = sp.expand(sp.series(G, t, 0, K + 2).removeO())
loes = {}
for k in range(1, K + 1):
    loes[cs[k - 1]] = sp.solve(sp.expand(G.coeff(t, k + 1).subs(loes)), cs[k - 1])[0]
c_fatou = [mp.mpf(sp.Rational(loes[cs[k - 1]]).p) / mp.mpf(sp.Rational(loes[cs[k - 1]]).q) for k in range(1, K + 1)]


def Phi(Z):
    return Z - mp.log(abs(Z)) / 3 + mp.fsum(c_fatou[k - 1] * Z ** (-k) for k in range(1, K + 1))


def Phi_strich(Z):
    return 1 - 1 / (3 * Z) - mp.fsum(k * c_fatou[k - 1] * Z ** (-k - 1) for k in range(1, K + 1))


def vor(y):
    return mp.exp(y / E)


def zurueck(y):
    return E * mp.log(y)


# Rückwärtsbahn von e² und Konstante
NMAX = 20000
bahn = [E ** 2]
for _ in range(NMAX):
    bahn.append(zurueck(bahn[-1]))
Cch = {n: Phi(2 * E / (E - bahn[n])) + n for n in (2500, 5000, 10000, 20000)}
lauf.pruefe("Konstante C_ch = lim Φ(Z(cheta(−n))) + n stabil (n = 2500 … 20000)", abs(Cch[2500] - Cch[20000]) < mp.mpf(10) ** -30,
            mp.nstr(Cch[20000], 30))
C_ch = Cch[20000]


def Z_aus(x):
    """Z < 0 mit Φ(Z) = x + C_ch (x stark negativ)."""
    Z = x + C_ch - mp.log(abs(x)) / 3
    for _ in range(60):
        dZ = (Phi(Z) - x - C_ch) / Phi_strich(Z)
        Z -= dZ
        if abs(dZ) < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
    return Z


def cheta(x, N=2000):
    y = E - 2 * E / Z_aus(mp.mpf(x) - N)
    for _ in range(N):
        y = vor(y)
    return y


ganz = {"−1": cheta(-1), "0": cheta(0), "1": cheta(1)}
lauf.pruefe("cheta(−1) = 2e, cheta(0) = e², cheta(1) = e^e (aus der Fatou-Koordinate zurückgewonnen)",
            abs(ganz["−1"] - 2 * E) + abs(ganz["0"] - E ** 2) + abs(ganz["1"] - E ** E) < mp.mpf(10) ** -25, {k: mp.nstr(v, 15) for k, v in ganz.items()})
h1, h2 = cheta(mp.mpf(1) / 2, N=1000), cheta(mp.mpf(1) / 2, N=3000)
lauf.pruefe("cheta(1/2) unabhängig von der Tiefe (N = 1000 gegen 3000)", abs(h1 - h2) < mp.mpf(10) ** -28, mp.nstr(h2, 25))
lauf.pruefe("Gleichung cheta(x+1) = η^cheta(x) bei x = −0.3", abs(cheta(mp.mpf("0.7")) - ETA ** cheta(mp.mpf("-0.3"))) < mp.mpf(10) ** -25)
asym = {n: (cheta(-n) - E) * n / (2 * E) for n in (100, 1000)}
lauf.pruefe("cheta(x) ≈ e + 2e/|x| für x → −∞", all(abs(v - 1) < 3 * mp.log(n) / n for n, v in asym.items()), {n: mp.nstr(v, 6) for n, v in asym.items()})

# untere η-Tetration zum Vergleich (C_T aus exp.eta_digamma)
C_T = mp.mpf("2.7982481542313876624525831854691003")


def T_eta(x, N=2000):
    Z = x + N + mp.log(x + N) / 3 + C_T
    for _ in range(60):
        dZ = (Phi(Z) - (x + N) - C_T) / Phi_strich(Z)
        Z -= dZ
        if abs(dZ) < mp.mpf(10) ** -40:
            break
    y = E - 2 * E / Z
    for _ in range(N):
        y = zurueck(y)
    return y


lauf.pruefe("Kontrolle untere Lösung: T_η(1/2) = 1.25715307505…", abs(T_eta(mp.mpf(1) / 2) - mp.mpf("1.25715307505417262217")) < mp.mpf(10) ** -18)

# ------------------------------------------------------------------ Szekeres: Abel-Funktion von e^u − 1 auf u > 0 ist total monoton
# Dieselbe Fatou-Reihe: Z = −2/u (u = y/e − 1). A(u) = Φ(−2/u_n) + n mit u_n = n-mal log(1 + ·) angewandt.
def A_sz(u, n=400):
    for _ in range(n):
        u = mp.log(1 + u)
    return Phi(-2 / u) + n


abl = {}
for x in (mp.mpf("0.05"), mp.mpf("0.3"), mp.mpf(1), mp.mpf(3), mp.mpf(10)):
    abl[mp.nstr(x, 3)] = [mp.diff(A_sz, x, k) for k in range(1, 9)]
lauf.pruefe("Szekeres: A′ > 0, A″ < 0, A‴ > 0, … bis zur 8. Ableitung (u = 0.05 … 10)",
            all((-1) ** (k + 1) * v[k - 1] > 0 for v in abl.values() for k in range(1, 9)),
            {x: [mp.nstr(v[k], 3) for k in (0, 1, 7)] for x, v in abl.items()})
lauf.pruefe("Abel-Gleichung A(e^u − 1) = A(u) + 1 (u = 0.7)", abs(A_sz(mp.expm1(mp.mpf("0.7"))) - A_sz(mp.mpf("0.7")) - 1) < mp.mpf(10) ** -30)
c_st = mp.mpf("0.001")


def A_gestoert(u):
    a = A_sz(u)
    return a + c_st * mp.sin(2 * mp.pi * a)


verletzt = None
for i in range(1, 200):
    x = mp.mpf(i) / 400
    if mp.diff(A_gestoert, x, 2) > 0:
        verletzt = x
        break
lauf.pruefe("gestörte Abel-Funktion A + 0.001·sin 2πA verletzt das schon bei der 2. Ableitung (nahe 0)", verletzt is not None,
            mp.nstr(verletzt, 3) if verletzt else None)

lauf.daten = {"usexp_sqrt2": {k: mp.nstr(mp.re(v), 30) for k, v in werte.items()}, "periode_sqrt2": mp.nstr(mp.im(per), 15),
              "C_ch": mp.nstr(C_ch, 30), "cheta_halb": mp.nstr(h2, 30), "cheta_minus_halb": mp.nstr(cheta(mp.mpf(-1) / 2), 30),
              "szekeres_ableitungen": {x: [mp.nstr(w, 6) for w in v] for x, v in abl.items()},
              "szekeres_stoerung_verletzt_bei": mp.nstr(verletzt, 4) if verletzt else None}

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
xs = [mp.mpf(i) / 10 for i in range(-60, 17)]
ax[0].plot([float(x) for x in xs if x > -1.9], [float(T_eta(x)) for x in xs if x > -1.9], color="#3d6fb6", label="untere Lösung T_η (T(0) = 1)")
ax[0].plot([float(x) for x in xs], [float(cheta(x)) for x in xs], color="#b5523b", label="obere Lösung cheta (cheta(0) = e²)")
ax[0].axhline(float(E), ls=":", color="#8a93a3"); ax[0].text(-5.8, float(E) + 0.2, "e", color="#8a93a3")
ax[0].set_ylim(-2, 16); ax[0].set_xlabel("x"); ax[0].set_title("Basis η: zwei Superexponentiale an e"); ax[0].legend(fontsize=8)
xs2 = [mp.mpf(i) / 10 for i in range(-80, 25)]
ax[1].plot([float(x) for x in xs2], [float(mp.re(usexp(x))) for x in xs2], color="#7a49c2", label="usexp_√2 (usexp(0) = 5)")
ax[1].axhline(4, ls=":", color="#8a93a3"); ax[1].text(-7.8, 4.15, "4", color="#8a93a3")
ax[1].set_ylim(3, 14); ax[1].set_xlabel("x"); ax[1].set_title("Basis √2: obere Superexponentialfunktion"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"usexp_√2(1/2) = {mp.nstr(mp.re(werte['1/2']), 15)} (ganz, Periode {mp.nstr(mp.im(per), 6)}i); "
                        f"cheta(1/2) = {mp.nstr(h2, 15)}, C_ch = {mp.nstr(C_ch, 15)}; cheta(x) − e ≈ 2e/|x| für x → −∞. "
                        f"Szekeres: Abel-Funktion von e^u − 1 total monoton bis Ordnung 8, Störung 0.001·sin 2πA verletzt das bei u = {mp.nstr(verletzt, 3)}.")
lauf.schreibe()
