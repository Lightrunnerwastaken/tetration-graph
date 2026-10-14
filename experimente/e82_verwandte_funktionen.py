"""exp.verwandte_funktionen — Produkt-Tetration (Superfunktion von x·aˣ, Forum tid 1808/1816) und exponentielle Fakultät aus sexp (tid 1828).

 1. Basiswechsel: g_a(x) = x·aˣ ist zu g_e konjugiert, ln a·g_a(x) = g_e(x·ln a); daher prodtet_a(c) = prodtet_e(c + prodslog_e(ln a))/ln a.
 2. 0 ist parabolischer Fixpunkt von x·e^{λx}; formale Iterierte b + cλb² + c(2c−1)/2·λ²b³ + …; x·eˣ und eˣ − 1 sind NICHT formal
    konjugiert (iterativer Rest ½ gegen ⅓) — die Produkt-Tetration ist eine echt andere parabolische Iteration.
 3. Natürliche (Fatou-)Produkt-Tetration zur Basis e auf der reellen Achse: prodtet(c) = α⁻¹(α(1) + c), α Fatou-Koordinate bei 0 (abstossende
    Seite), Abstieg mit Lambert W. Vergleich mit dem Interpolationswert 1.5134 aus dem Forum bei c = 1/2.
 4. Exponentielle Fakultät: D = lim (slog_e(n!₂) − n) ≈ −0.0613312423001 (Forum) und Fortsetzung
    x!₂ = lim_k log_{x+1}(log_{x+2}(… log_{x+k}(sexp(x+k+D)))) — Ganzzahlprobe und Konvergenz bei x = 1/2.
"""
import mpmath as mp
import sympy as sp

from gemeinsam import Lauf, cache_speichern, kneser, plt

mp.mp.dps = 40
lauf = Lauf("exp.verwandte_funktionen", __file__)
daten = {}
X = sp.Symbol("x")

# ------------------------------------------------------------------ 1. Basiswechsel
d = max(abs(mp.log(a) * x * a ** x - (x * mp.log(a)) * mp.exp(x * mp.log(a))) for a in (mp.mpf(2), mp.mpf(10)) for x in (mp.mpf("0.3"), mp.mpf(2)))
lauf.pruefe("ln a · x·aˣ = (x ln a)·e^{x ln a}: x·aˣ ist zu x·eˣ konjugiert (Streckung um ln a)", d < mp.mpf(10) ** -35, mp.nstr(d, 3))

# ------------------------------------------------------------------ 2. formale Iteration am parabolischen Punkt 0
lam, c = sp.symbols("lambda c")
f = sp.series(X * sp.exp(lam * X), X, 0, 5).removeO()
h = X + sp.Rational(1, 2) * lam * X ** 2 + 0 * X ** 3          # c = 1/2: Koeffizient c(2c−1)/2·λ² = 0
hh = sp.expand(h.subs(X, h))
diff3 = sp.expand(sp.series(hh - f, X, 0, 4).removeO())
lauf.pruefe("Halbe Iterierte von x·e^{λx}: h = x + λx²/2 + 0·x³ + …, h∘h = x·e^{λx} bis x³ (Formel c(2c−1)/2·λ² bei c = 1/2)", diff3 == 0, str(diff3))

ORD = 6
F1 = sp.series(X * sp.exp(X), X, 0, ORD + 1).removeO()
F2 = sp.series(sp.exp(X) - 1, X, 0, ORD + 1).removeO()
rest = lambda F: 1 - F.coeff(X, 3) / F.coeff(X, 2) ** 2   # iterativer Rest von x + a₂x² + a₃x³ + … (Formel-Invariante)
lauf.pruefe("x·eˣ und eˣ − 1 sind NICHT formal konjugiert (auch nicht mit Streckung): die Formel-Invariante 1 − a₃/a₂² ist ½ gegen ⅓",
            rest(F1) == sp.Rational(1, 2) and rest(F2) == sp.Rational(1, 3), {"x·eˣ": str(rest(F1)), "eˣ−1": str(rest(F2))})

# ------------------------------------------------------------------ 3. Fatou-Koordinate von x·eˣ und Produkt-Tetration
K = 6
ds = sp.symbols(f"d1:{K + 1}")
alpha_r = -1 / X + sp.Rational(1, 2) * sp.log(X) + sum(ds[i] * X ** (i + 1) for i in range(K))
fx = sp.series(X * sp.exp(X), X, 0, K + 4).removeO()
ausdruck = sp.series(-1 / fx + sp.Rational(1, 2) * sp.log(sp.expand(fx / X)) + sum(ds[i] * fx ** (i + 1) for i in range(K))
                     - (-1 / X + sum(ds[i] * X ** (i + 1) for i in range(K))) - 1, X, 0, K + 1).removeO()
ausdruck = sp.expand(ausdruck)
dl = sp.solve([ausdruck.coeff(X, k) for k in range(1, K + 1)], ds, dict=True)[0]
lauf.pruefe("Fatou-Koordinate α(x) = −1/x + ½ ln x + Σ d_k x^k löst α(x eˣ) = α(x) + 1 formal (½ = 1 − a₃/a₂², a₃ = ½)",
            sp.expand(ausdruck.coeff(X, 0)) == 0, {str(k): str(v) for k, v in list(dl.items())[:3]})
afun = sp.lambdify(X, alpha_r.subs(dl).subs({d: 0 for d in ds}), "mpmath")
g = lambda t: t * mp.exp(t)
W = lambda t: mp.lambertw(t).real
NAB = 20000


T0 = mp.mpf(1)
for _ in range(NAB):
    T0 = W(T0)                                   # W^{NAB}(1), klein; α(1) = α_asym(T0) + NAB


def prodtet(cc):
    """α(prodtet(c)) = α(1) + c: in der Nähe von T0 in der Fatou-Koordinate um c verschieben, dann NAB-mal vorwärts."""
    ziel = afun(T0) + cc
    t0 = T0 + cc * T0 ** 2                       # α′ ≈ 1/t², also Δt ≈ c·t²
    t = mp.findroot(lambda s: afun(s) - ziel, (t0, t0 * (1 + mp.mpf(10) ** -8)))
    for _ in range(NAB):
        t = g(t)
    return t


p = {str(cc): prodtet(mp.mpf(cc)) for cc in ("-1", "0", "0.5", "1")}
lauf.pruefe("Fatou-Produkt-Tetration: prodtet(0) = 1, prodtet(1) = e, prodtet(−1) = W(1) = 0.567143…",
            abs(p["0"] - 1) < 1e-20 and abs(p["1"] - mp.e) < 1e-19 and abs(p["-1"] - W(1)) < 1e-20,
            {k: mp.nstr(v, 15) for k, v in p.items()})
halbe = p["0.5"]
daten["prodtet_e_halb"] = mp.nstr(halbe, 20)
lauf.pruefe("prodtet_e(½) aus der Fatou-Koordinate = 1.513428… — stimmt mit dem Interpolationswert 1.5134 aus dem Forum auf allen angegebenen Stellen",
            abs(halbe - mp.mpf("1.5134")) < mp.mpf("0.00005"), {"Fatou": mp.nstr(halbe, 12), "Forum (Lagrange)": "1.5134"})
cs_ = [mp.mpf(j) / 10 - 1 for j in range(31)]
kurve = [prodtet(cc) for cc in cs_]

# ------------------------------------------------------------------ 4. exponentielle Fakultät und D
def sexp_e(t):
    """Kneser-sexp zur Basis e; für grosse t über exp aus dem Grundintervall."""
    t = mp.mpf(t)
    m = int(mp.floor(t + 1)) if t > 0 else 0
    v = mp.re(kneser("sexp", "exp(1)", [t - m], dps=40)[0])
    for _ in range(m):
        v = mp.exp(v)
    return v


def slog_e(Y):
    Y, k = mp.mpf(Y), 0
    while Y > 2:
        Y, k = mp.log(Y), k + 1
    return k + mp.re(kneser("slog", "exp(1)", [Y], dps=40)[0])


EF = {1: mp.mpf(1)}
for n in range(2, 6):
    EF[n] = mp.mpf(n) ** EF[n - 1]                  # EF(5) = 5^262144 ist als mpf darstellbar
Dn = {n: slog_e(EF[n]) - n for n in range(2, 6)}
Dn[6] = slog_e(EF[5] * mp.log(6)) + 1 - 6           # ln EF(6) = EF(5)·ln 6
Dn[7] = slog_e(EF[5] * mp.log(6) + mp.log(mp.log(7))) + 2 - 7   # ln ln EF(7) = ln EF(6) + ln ln 7
daten["D_n"] = {n: mp.nstr(v, 15) for n, v in Dn.items()}
D = Dn[7]
lauf.pruefe("D = lim slog_e(n!₂) − n = −1.0552140087 (Kneser, slog(1) = 0), stabil ab n = 6",
            abs(Dn[7] - Dn[6]) < mp.mpf("1e-13") and abs(Dn[7] - mp.mpf("-1.0552140087127")) < mp.mpf("1e-12"),
            {n: mp.nstr(v, 14) for n, v in Dn.items()})


def ef_fort(x, k, DD=None):
    """log_{x+1}(log_{x+2}(… log_{x+k}(sexp(x+k+D)))); innen eine Stufe logarithmisch: ln log_{x+k}(sexp(t)) = sexp(t−2) − ln ln(x+k)."""
    x = mp.mpf(x)
    L = sexp_e(x + k + (D if DD is None else DD) - 2) - mp.log(mp.log(x + k))
    for j in range(k - 1, 0, -1):
        v = L / mp.log(x + j)
        if j == 1:
            return v
        L = mp.log(v)


tab = {"1": [ef_fort("1", k) for k in (5, 6)], "2": [ef_fort("2", k) for k in (4, 5)],
       "0.5": [ef_fort("0.5", k) for k in (5, 6)], "1.5": [ef_fort("1.5", k) for k in (4, 5)]}
cache_speichern()
daten["fortsetzung"] = {x0: [mp.nstr(v, 15) for v in vs] for x0, vs in tab.items()}
lauf.pruefe("Ganzzahlprobe: die Fortsetzung gibt 1!₂ = 1 und 2!₂ = 2 auf 10⁻¹²",
            abs(tab["1"][-1] - 1) < 1e-12 and abs(tab["2"][-1] - 2) < 1e-12, {"1": mp.nstr(tab["1"][-1], 15), "2": mp.nstr(tab["2"][-1], 15)})
lauf.pruefe("bei x = ½ und 1.5 konvergiert die Kette (letzte Verlängerung ändert < 10⁻¹²)",
            abs(tab["0.5"][-1] - tab["0.5"][-2]) < 1e-12 and abs(tab["1.5"][-1] - tab["1.5"][-2]) < 1e-12,
            {"0.5!₂": mp.nstr(tab["0.5"][-1], 15), "1.5!₂": mp.nstr(tab["1.5"][-1], 15)})
forum = ef_fort("1", 6, mp.mpf("-0.0613312423001"))
daten["mit_forum_D"] = mp.nstr(forum, 12)
lauf.pruefe("Der Forumswert D = −0.0613312423001 passt nicht zu Kneser mit slog(1) = 0: die Formel gäbe damit 1!₂ ≈ 3.5",
            abs(forum - 1) > 1, mp.nstr(forum, 12))

nahe0 = {x0: ef_fort(x0, 7) for x0 in ("0.1", "0.01", "0.001", "0.0001", "0.00001")}
daten["gegen_0"] = {x0: mp.nstr(v, 15) for x0, v in nahe0.items()}
steigung = (nahe0["0.0001"] - nahe0["0.00001"]) / mp.mpf("0.00009")
grenz0 = nahe0["0.00001"] - steigung * mp.mpf("0.00001")
daten["EF_0"] = mp.nstr(grenz0, 10)
lauf.pruefe("TPID 2 an dieser Fortsetzung: x!₂ → 0.57305 für x ↓ 0 — nahe an γ = 0.57722, aber um 0.004 daneben",
            abs(nahe0["0.0001"] - nahe0["0.00001"]) < mp.mpf("1e-4") and mp.mpf("0.003") < abs(grenz0 - mp.euler) < mp.mpf("0.006"),
            {x0: mp.nstr(v, 10) for x0, v in nahe0.items()} | {"Grenzwert": mp.nstr(grenz0, 8), "gamma": mp.nstr(mp.euler, 10)})

lauf.daten = daten
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
a1.plot([float(t) for t in cs_], [float(v) for v in kurve])
a1.set_yscale("log"); a1.set_xlabel("c"); a1.set_title("Produkt-Tetration zur Basis e (Fatou-Koordinate)")
xs = [mp.mpf(j) / 20 for j in range(1, 46)]
a2.plot([float(t) for t in xs], [float(ef_fort(t, 5 if t < 1.6 else 4)) for t in xs], label="x!₂ über sexp")
a2.axhline(float(mp.euler), ls=":", c="gray", label="γ")
a2.plot([1, 2], [1, 2], "o", label="ganze Werte")
a2.set_xlabel("x"); a2.set_title("Exponentielle Fakultät aus sexp"); a2.legend(fontsize=8)
lauf.bild(fig)
cache_speichern()
lauf.zusammenfassung = (f"Produkt-Tetration: Basiswechsel als Konjugation, parabolischer Fixpunkt 0 (nicht formal konjugiert zu eˣ − 1), "
                        f"Fatou-Wert prodtet_e(½) = {mp.nstr(halbe, 8)} (Forum 1.5134). Exponentielle Fakultät: D = {mp.nstr(D, 13)} "
                        f"(Forumswert −0.0613 nicht reproduziert), Fortsetzung konvergiert, 0.5!₂ = {mp.nstr(tab['0.5'][-1], 10)}, "
                        f"x!₂ → {mp.nstr(grenz0, 6)} für x ↓ 0 (γ = 0.57722).")
lauf.schreibe()
