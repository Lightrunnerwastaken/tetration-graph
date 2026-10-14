"""exp.stoerungsreihe — Tetration als Störung von eˣ − 1 (Polynom ohne Carleman, Forum tid 1827) und eine Näherungsleiter für Desmos.

 1. Einbettung: G(y) = e^{G(y−1)} − 1 + u, G(0) = 0. Bei u = 1 ist G(y) = sexp_e(y − 1); bei u = 0 bleibt G ≡ 0 (parabolischer Fixpunkt von eˣ − 1).
    Entwicklung G = Σ g_n(y) uⁿ: g_1 = y, g_n(y) − g_n(y−1) = [uⁿ](e^{G(y−1)} − 1 − G(y−1)) hängt nur von g_1 … g_{n−1} ab, also ist g_n ein
    Polynom vom Grad 2n − 1 (Faulhaber-Summe). Das ist die Rekursion des Forums mit f_n = n!·g_n (Bell-Polynome = Koeffizienten von e^G).
 2. P_N(x) = Σ_{n≤N} g_n(x+1): N = 8 gibt das Polynom vom Grad 15 aus dem Forum (a₀ = 1, a₁ = 87159883/79833600).
 3. Abstand zu Kneser auf [−1, 0] in Abhängigkeit von N: die Reihe in u divergiert, es gibt eine beste Ordnung.
 4. Näherungsleiter für Desmos: linear, P_N, Tschebyschow-Interpolation an Kneser; Fehler und Knick (Sprung von T′ bei 0) der stückweisen Fortsetzung.
"""
import mpmath as mp
import sympy as sp

from gemeinsam import Lauf, cache_speichern, gauss_legendre, kneser, plt

mp.mp.dps = 40
lauf = Lauf("exp.stoerungsreihe", __file__)
daten = {}
y = sp.Symbol("y")
NMAX = 14


QQ = sp.QQ
BERN = {}


def faulhaber(q):
    """S mit S(y) − S(y−1) = q(y−1), S(0) = 0, also S(y) = Σ_{t=0}^{y−1} q(t)."""
    s = sp.Poly(0, y, domain=QQ)
    for (m,), c in q.terms():
        if m not in BERN:
            BERN[m] = sp.Poly((sp.bernoulli(m + 1, y) - sp.bernoulli(m + 1, 0)) / (m + 1), y, domain=QQ)
        s += BERN[m] * c
    return s


def verschoben(p, a):
    return sp.Poly(p.as_expr().subs(y, y + a), y, domain=QQ)


def exp_koef(gs, n):
    """Koeffizienten E_0 … E_n von exp(Σ gs[k] u^k) (gs[0] unbenutzt), über n E_n = Σ k g_k E_{n−k}."""
    E = [sp.Poly(1, y, domain=QQ)]
    for k in range(1, n + 1):
        E.append(sum((gs[j] * E[k - j] * j for j in range(1, k + 1) if j < len(gs)), sp.Poly(0, y, domain=QQ)) * QQ(1, k))
    return E


g = [None, sp.Poly(y, y, domain=QQ)]            # g[n](y)
for n in range(2, NMAX + 1):
    gm = [None] + [verschoben(gk, -1) for gk in g[1:]]                  # g_k(y − 1), nur k < n
    q = verschoben(exp_koef(gm, n)[n], 1)        # [uⁿ] e^{G(y−1)} ohne g_n-Anteil, als q(y−1)
    g.append(faulhaber(q))

grade = [g[n].degree() for n in range(1, NMAX + 1)]
lauf.pruefe("g_n ist ein Polynom vom Grad 2n − 1 (n = 1 … 14)", grade == [2 * n - 1 for n in range(1, NMAX + 1)], grade)

# Funktionalgleichung an gebrochenen Höhen: [uⁿ] G(y) = [uⁿ] (e^{G(y−1)} − 1 + u), exakt in ℚ
ok = True
for yv in (sp.Rational(1, 2), sp.Rational(7, 3)):
    links = [None] + [gk.eval(yv) for gk in g[1:9]]
    unten = [None] + [sp.Poly(gk.eval(yv - 1), y, domain=QQ) for gk in g[1:9]]
    E = exp_koef(unten, 8)
    rechts = [None, E[1].eval(0) + 1] + [E[k].eval(0) for k in range(2, 9)]
    ok &= all(links[k] == rechts[k] for k in range(1, 9))
lauf.pruefe("Funktionalgleichung G(y) = e^{G(y−1)} − 1 + u gilt exakt bis u⁸ auch bei y = 1/2 und 7/3", ok)
g = [None] + [gk.as_expr() for gk in g[1:]]

x = sp.Symbol("x")
P8 = sp.Poly(sp.expand(sum(g[n].subs(y, x + 1) for n in range(1, 9))), x)
koef = P8.all_coeffs()[::-1]
daten["P8_koeffizienten"] = [str(c) for c in koef]
lauf.pruefe("N = 8 gibt das Forumspolynom: Grad 15, a₀ = 1, a₁ = 87159883/79833600",
            P8.degree() == 15 and koef[0] == 1 and koef[1] == sp.Rational(87159883, 79833600), [str(c) for c in koef[:3]])


def P(N):
    f = sp.lambdify(x, sp.expand(sum(g[n].subs(y, x + 1) for n in range(1, N + 1))), "mpmath")
    return lambda t: f(mp.mpf(t))


def Pd(N):
    f = sp.lambdify(x, sp.diff(sum(g[n].subs(y, x + 1) for n in range(1, N + 1)), x), "mpmath")
    return lambda t: f(mp.mpf(t))


ks = lambda zs: [mp.re(v) for v in kneser("sexp", "exp(1)", zs, dps=40)]

# ------------------------------------------------------------------ 3. Abstand zu Kneser in Abhängigkeit von N
fehler = {}
for N in range(1, NMAX + 1):
    pn = P(N)
    fehler[N] = gauss_legendre(lambda zs: [abs(k - pn(z)) for z, k in zip(zs, ks(zs))], mp.mpf(-1), mp.mpf(0), 40)
daten["L1_fehler_nach_N"] = {N: mp.nstr(v, 4) for N, v in fehler.items()}
Nbest = min(fehler, key=fehler.get)
trapez = None                                            # Gegenprobe mit 2001-Punkte-Trapez weiter unten
lauf.pruefe("Die Reihe divergiert: der Fehler fällt bis zu einer besten Ordnung und steigt danach wieder",
            fehler[NMAX] > fehler[Nbest] * 10 and 6 <= Nbest <= 10, {"beste_N": Nbest, "Fehler_bestes": mp.nstr(fehler[Nbest], 3), f"N={NMAX}": mp.nstr(fehler[NMAX], 3)})
mitte = {N: mp.nstr(g[N].subs(y, sp.Rational(1, 2)), 5) for N in range(1, NMAX + 1)}
daten["g_n(1/2)"] = mitte

# ------------------------------------------------------------------ 4. Näherungsleiter (Desmos)
gitter = [mp.mpf(-1) + mp.mpf(j) / 2000 for j in range(2001)]
kg = ks(gitter)
p8 = P(8)
e8 = [abs(k - p8(t)) for t, k in zip(gitter, kg)]
trapez = (mp.fsum(e8) - (e8[0] + e8[-1]) / 2) / 2000
lauf.pruefe("∫₋₁⁰ |sexp − P_8| = 2.4·10⁻⁷ mit zwei Quadraturen (Forumsangabe 5.52·10⁻⁷: gleiche Grössenordnung, aber doppelt so gross)",
            abs(trapez - fehler[8]) < mp.mpf("0.05e-7") and mp.mpf("1e-7") < fehler[8] < mp.mpf("5.52e-7"),
            {"Gauss-Legendre 40": mp.nstr(fehler[8], 4), "Trapez 2001": mp.nstr(trapez, 4)})


def leiter_eintrag(name, f, fd, grad):
    maxf = max(abs(k - f(t)) for t, k in zip(gitter, kg))
    knick = abs(fd(mp.mpf(0)) - fd(mp.mpf(-1)))           # T′(0⁻) = f′(0), T′(0⁺) = e^{f(−1)} f′(−1) = f′(−1), da f(−1) = 0
    f2 = lambda t: mp.diff(f, t, 2)
    knick2 = abs(f2(mp.mpf(0)) - (f2(mp.mpf(-1)) + fd(mp.mpf(-1)) ** 2))   # T″(0⁺) = f″(−1) + f′(−1)²
    return {"name": name, "grad": grad, "max_fehler": mp.nstr(maxf, 3), "knick_bei_0": mp.nstr(knick, 3),
            "sprung_T2_bei_0": mp.nstr(knick2, 3), "_f": float(maxf), "_k": float(knick)}


leiter = [leiter_eintrag("linear 1 + x", lambda t: 1 + t, lambda t: mp.mpf(1), 1)]
for N in (2, 4, 8):
    leiter.append(leiter_eintrag(f"Störungsreihe P_{N}", P(N), Pd(N), 2 * N - 1))


def tscheb(d):
    """Interpolation an Kneser in d + 1 Tschebyschow-Knoten von [−1, 0], dazu die Randbedingungen f(−1) = 0, f(0) = 1 exakt."""
    knoten = [mp.mpf(-1)] + [(-1 + mp.cos(mp.pi * (2 * j + 1) / (2 * (d - 1)))) / 2 for j in range(d - 1)] + [mp.mpf(0)]
    werte = [mp.mpf(0)] + ks(knoten[1:-1]) + [mp.mpf(1)]
    A = mp.matrix([[t ** i for i in range(d + 1)] for t in knoten])
    c = mp.lu_solve(A, mp.matrix(werte))
    return (lambda t: mp.fsum(c[i] * t ** i for i in range(d + 1)),
            lambda t: mp.fsum(i * c[i] * t ** (i - 1) for i in range(1, d + 1)), [mp.nstr(c[i], 12) for i in range(d + 1)])


tsch = {}
for d in (3, 5, 8, 12):
    f, fd, c = tscheb(d)
    tsch[d] = c
    leiter.append(leiter_eintrag(f"Tschebyschow an Kneser, Grad {d}", f, fd, d))
cache_speichern()
daten["leiter"] = [{k: v for k, v in e.items() if not k.startswith("_")} for e in leiter]
lauf.pruefe("Die Störungsreihe ist an den Nahtstellen exakt C¹ (jedes g_n hat g_n′(1) = g_n′(0)), aber nicht C²",
            all(e["_k"] < 1e-30 for e in leiter[1:4]), {e["name"]: e["sprung_T2_bei_0"] for e in leiter[1:4]})
daten["tschebyschow_koeffizienten"] = tsch
fl = {e["name"]: e["_f"] for e in leiter}
lauf.pruefe("Leiter: linear ~10⁻², P_8 ~10⁻⁶, Tschebyschow Grad 8 < 10⁻⁶, Grad 12 < 10⁻⁸",
            1e-3 < fl["linear 1 + x"] < 1e-1 and fl["Störungsreihe P_8"] < 1e-5
            and fl["Tschebyschow an Kneser, Grad 8"] < 1e-6 and fl["Tschebyschow an Kneser, Grad 12"] < 1e-8,
            {k: f"{v:.2e}" for k, v in fl.items()})

lauf.daten = daten
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
a1.semilogy(list(fehler), [float(v) for v in fehler.values()], "o-")
a1.set_xlabel("Ordnung N"); a1.set_ylabel("∫₋₁⁰ |sexp − P_N|"); a1.set_title("Störungsreihe: beste Ordnung, dann Divergenz")
a2.semilogy([e["grad"] for e in leiter[1:4]], [e["_f"] for e in leiter[1:4]], "o-", label="Störungsreihe P_N")
a2.semilogy([e["grad"] for e in leiter[4:]], [e["_f"] for e in leiter[4:]], "s-", label="Tschebyschow an Kneser")
a2.axhline(leiter[0]["_f"], ls="--", c="gray", label="linear 1 + x")
a2.set_xlabel("Polynomgrad"); a2.set_ylabel("max. Fehler auf [−1, 0]"); a2.set_title("Näherungsleiter"); a2.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Forumspolynom (tid 1827) als Störungsreihe von eˣ − 1 hergeleitet und exakt nachgebaut; ∫|Fehler| bei N = 8 = "
                        f"{mp.nstr(fehler[8], 3)}, beste Ordnung N = {Nbest}, danach Divergenz. Näherungsleiter für Desmos gemessen.")
lauf.schreibe()
