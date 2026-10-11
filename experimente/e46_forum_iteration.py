"""exp.forum_iteration — Iterations-Themen aus dem Forum: elementare Superfunktionen, die Möbius-Familie über den parabolischen Fall,
 TPID 9 (intuitiver Logarithmus, Trappmann 2010) und TPID 13 (Newton-Interpolation von n^(1/n)).
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 60
lauf = Lauf("exp.forum_iteration", __file__)

# 1. Elementare Superfunktionen F(t+1) = f(F(t)) an zufälligen Stellen
c = mp.mpf("0.37")
paare = {
    "x² − 2 ↔ 2cosh(2^t c)": (lambda x: x ** 2 - 2, lambda t: 2 * mp.cosh(mp.power(2, t) * c)),
    "2x² − 1 ↔ cosh(2^t c)": (lambda x: 2 * x ** 2 - 1, lambda t: mp.cosh(mp.power(2, t) * c)),
    "(x−1)² + 1 ↔ 2^(2^t) + 1 (Fermat)": (lambda x: (x - 1) ** 2 + 1, lambda t: mp.power(2, mp.power(2, t)) + 1),
    "x³ + 3x ↔ 2sinh(3^t c)": (lambda x: x ** 3 + 3 * x, lambda t: 2 * mp.sinh(mp.power(3, t) * c)),
    "4x(1 − x) ↔ sin²(2^t c)": (lambda x: 4 * x * (1 - x), lambda t: mp.sin(mp.power(2, t) * c) ** 2),
    "x² + 2x ↔ e^(2^t c) − 1": (lambda x: x ** 2 + 2 * x, lambda t: mp.exp(mp.power(2, t) * c) - 1),
    "x^a ↔ e^(a^t c) (a = 1.7)": (lambda x: mp.power(x, mp.mpf("1.7")), lambda t: mp.exp(mp.power(mp.mpf("1.7"), t) * c)),
    "2x/(1+x) ↔ 2^t/(2^t + c)": (lambda x: 2 * x / (1 + x), lambda t: mp.power(2, t) / (mp.power(2, t) + c)),
    "x/(1+x) ↔ 1/(t + c)": (lambda x: x / (1 + x), lambda t: 1 / (t + c)),
}
res = 0
for name, (f, F) in paare.items():
    for t in (mp.mpf("-0.3"), mp.mpf("0.41"), mp.mpf("1.7")):
        res = max(res, abs(F(t + 1) - f(F(t))) / abs(F(t + 1)))
lauf.pruefe("neun elementare Superfunktionen erfüllen F(t+1) = f(F(t)) (relativ)", res < mp.mpf(10) ** -50, mp.nstr(res, 3))

# 2. Möbius-Familie A(ε) = [[1, 1], [−ε, 1]]: f^[t] über A^t ist in ε analytisch, auch durch den parabolischen Fall ε = 0
def halb(eps, z, t=mp.mpf("0.5")):
    A = mp.matrix([[1, 1], [-eps, 1]])
    At = mp.expm(t * mp.logm(A))
    return (At[0, 0] * z + At[0, 1]) / (At[1, 0] * z + At[1, 1])


def f(eps, z):
    return (z + 1) / (-eps * z + 1)


z0 = mp.mpf("0.3")
werte = {e: halb(mp.mpf(e), z0) for e in ("-0.01", "-0.001", "0", "0.001", "0.01")}
ok_halb = all(abs(halb(mp.mpf(e), halb(mp.mpf(e), z0)) - f(mp.mpf(e), z0)) < mp.mpf(10) ** -45 for e in werte)
# Analytizität in ε: zweite Differenz glatt (symmetrische Ableitungen von beiden Seiten gleich)
h = mp.mpf("0.001")
d_links = (werte["0"] - werte["-0.001"]) / h
d_rechts = (werte["0.001"] - werte["0"]) / h
d_glatt = mp.diff(lambda e: mp.re(halb(e, z0)), 0)
lauf.pruefe("Möbius-Familie: halbe Iterierte reell für ε < 0, = 0, > 0 und glatt in ε durch den parabolischen Fall",
            ok_halb and all(abs(mp.im(w)) < mp.mpf(10) ** -40 for w in werte.values()) and abs((d_links + d_rechts) / 2 - d_glatt) < 1e-5,
            {e: mp.nstr(mp.re(w), 12) for e, w in werte.items()})

# 3. TPID 9: a'_1 = 1, a'_n = 1/(1 − b^n) Σ_{m<n} a'_m C(n,m)(1−b)^(n−m) b^m → (b−1)/ln b
def tpid9(b, N):
    a = [None, mp.mpf(1)]
    for n in range(2, N + 1):
        a.append(mp.fsum(a[m] * mp.binomial(n, m) * (1 - b) ** (n - m) * b ** m for m in range(1, n)) / (1 - b ** n))
    return a


b = mp.mpf(2)
a = tpid9(b, 120)
ziel = (b - 1) / mp.log(b)
lauf.pruefe("TPID 9 (Trappmann 2010): a'_n → (b−1)/ln b für b = 2, langsam und schwankend (|Abstand| < 5·10⁻⁵ ab n = 80)", all(abs(a[n] - ziel) < mp.mpf("5e-5") for n in range(80, 121)),
            {n: mp.nstr(a[n] - ziel, 3) for n in (10, 40, 80, 120)})

# 4. TPID 13: Newton-Interpolation von (0,0), (m, m^(1/m))
werte13 = [mp.mpf(0)] + [mp.power(m, mp.mpf(1) / m) for m in range(1, 90)]


def newton(a, x, N):
    s, tab = mp.mpf(0), {}
    for n in range(N + 1):
        s += mp.binomial(x, n) * mp.fsum((-1) ** (n - j) * mp.binomial(n, j) * a[j] for j in range(n + 1))
        if n in (20, 40, 80):
            tab[n] = s
    return tab


mp.mp.dps = 120
tab = {}
for x in (mp.mpf("2.5"), mp.mpf("7.5")):
    t = newton(werte13, x, 80)
    tab[mp.nstr(x, 2)] = {n: mp.nstr(abs(v - mp.power(x, 1 / x)), 3) for n, v in t.items()}
lauf.pruefe("TPID 13: Newton-Interpolation von m^(1/m) nähert x^(1/x) bei x = 2.5 und 7.5 (Abstand fällt mit N)",
            all(mp.mpf(v[80]) < mp.mpf(v[20]) for v in tab.values()), tab)
mp.mp.dps = 60

lauf.daten = {"moebius_halb": {e: mp.nstr(mp.re(w), 15) for e, w in werte.items()}, "tpid13": tab}
fig, ax = plt.subplots(1, 1, figsize=(6, 3.6))
es = [mp.mpf(i) / 200 for i in range(-40, 41)]
ax.plot([float(e) for e in es], [float(mp.re(halb(e, z0))) for e in es])
ax.axvline(0, ls=":", color="#888")
ax.set_xlabel("ε (ε < 0: zwei reelle Fixpunkte, ε > 0: komplexes Paar)"); ax.set_title("Möbius-Modell: halbe Iterierte glatt durch den parabolischen Fall")
lauf.bild(fig)
lauf.zusammenfassung = "Neun elementare Superfunktionen bestätigt; Möbius-Familie analytisch durch den parabolischen Fall; TPID 9 konvergiert; TPID 13 numerisch plausibel."
lauf.schreibe()
