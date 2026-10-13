"""exp.superwurzel_omega — Formeln aus „On Tetration I/II“ (What the Hectogon?, 2020) nachgerechnet.

 1. ε_n(x) = x·(eˣ)↑↑n (ε_0 = x, ε_{n+1} = x·e^{ε_n}), ω_n = Umkehrung: die n-te Superwurzel ist exp(ω_{n−1}(ln Y)).
 2. Unendlicher Turm: h(x) = 1/ssrt(1/x); ε_∞(x) = −W(−x), ω_∞(x) = x·e^{−x}.
 3. Basisableitung ohne Rekursion: (x↑↑n)' = x↑↑n/(x ln x) · Σ_{k=1}^n Π_{m=k−1}^{n−1} (x↑↑m · ln x).
 4. Ableitungen: ε_n' = (ε_n/x)·(1 + S_n(x)), ω_n' = ω_n/(x·(1 + S_n(ω_n))), S_n(x) = Σ_{k=1}^n Π_{m=k−1}^{n−1} ε_m(x).
 5. Grenzidentität: S_n(ω_n(x)) → x/(1−x) für 0 < x < 1 — aber NICHT, wenn man vorher ω_∞ einsetzt (Video): Σ_k Π_{m≥k−1} ε_m(x e^{−x}) → x.
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 40
lauf = Lauf("exp.superwurzel_omega", __file__)
daten = {}
TOL = mp.mpf(10) ** -30


def T(n, x):
    t = mp.mpf(0) if n == -1 else mp.mpf(1)
    for _ in range(max(n, 0)):
        t = x ** t
    return t


def eps(n, x):
    e = x
    for _ in range(n):
        e = x * mp.exp(e)
    return e


def omega(n, y, start=mp.mpf("0.3")):
    return mp.findroot(lambda x: eps(n, x) - y, start)


def omega_klein(n, y):                    # für 0 < y < 1: Wurzel liegt in (0, 1/e), dort konvergiert ε_n
    return mp.findroot(lambda x: eps(n, x) - y, (mp.mpf(0), mp.mpf("0.3678")), solver="anderson")


def S(n, x):
    return mp.fsum(mp.fprod(eps(m, x) for m in range(k - 1, n)) for k in range(1, n + 1))


# ------------------------------------------------------------------ 1. Superwurzeln
werte = {}
fehler = []
for n, Y in ((2, 2), (3, 2), (3, 10), (4, 3), (5, 10)):
    y = mp.exp(omega(n - 1, mp.log(Y), mp.mpf("0.5")))
    werte[f"{n}. Superwurzel von {Y}"] = mp.nstr(y, 25)
    fehler.append(abs(T(n, y) - Y))
w2 = mp.exp(mp.lambertw(mp.log(2)).real)
lauf.pruefe("n-te Superwurzel = exp(ω_{n−1}(ln Y)); n = 2 gibt e^{W(ln 2)}; dritte Superwurzel von 2 = 1.4766843373578699…",
            max(fehler) < TOL and abs(mp.mpf(werte["2. Superwurzel von 2"]) - w2) < mp.mpf(10) ** -23
            and werte["3. Superwurzel von 2"].startswith("1.47668433735786994708"), werte)
daten["superwurzeln"] = werte

# ------------------------------------------------------------------ 2. unendlicher Turm, ε_∞, ω_∞
d = []
for x in (mp.mpf("1.1"), mp.mpf("1.3"), mp.mpf("1.44")):
    h = mp.lambertw(-mp.log(x)).real / -mp.log(x)
    ssrt = lambda Y: mp.log(Y) / mp.lambertw(mp.log(Y)).real
    d.append(abs(h - 1 / ssrt(1 / x)))
y = mp.mpf("0.2")
e_inf = eps(400, y)
lauf.pruefe("h(x) = 1/ssrt(1/x) (x = 1.1, 1.3, 1.44); ε_∞(0.2) = −W(−0.2); ω_∞(x) = x e^{−x} kehrt ε_∞ um",
            max(d) < TOL and abs(e_inf + mp.lambertw(-y).real) < TOL
            and abs(-mp.lambertw(-(y * mp.exp(-y))).real - y) < TOL,
            {"max |h − 1/ssrt(1/x)|": mp.nstr(max(d), 3), "ε_400(0.2)": mp.nstr(e_inf, 20)})

# ------------------------------------------------------------------ 3. Basisableitung geschlossen
d = []
for n in (2, 3, 4, 5):
    for x in (mp.mpf("1.2"), mp.mpf("1.4")):
        f = T(n, x) / (x * mp.log(x)) * mp.fsum(mp.fprod(T(m, x) * mp.log(x) for m in range(k - 1, n)) for k in range(1, n + 1))
        d.append(abs(f - mp.diff(lambda t: T(n, t), x)))
lauf.pruefe("(x↑↑n)' = x↑↑n/(x ln x)·Σ_k Π_{m=k−1}^{n−1}(x↑↑m·ln x) für n = 2…5, x = 1.2, 1.4", max(d) < mp.mpf(10) ** -25,
            mp.nstr(max(d), 3))

# ------------------------------------------------------------------ 4. Ableitungen von ε_n und ω_n
d = []
for n in (1, 2, 3, 4):
    x = mp.mpf("0.3")
    d.append(abs(mp.diff(lambda t: eps(n, t), x) - eps(n, x) / x * (1 + S(n, x))))
    yv = mp.mpf("0.8")
    w = omega(n, yv)
    d.append(abs(mp.diff(lambda t: omega(n, t, w), yv) - w / (yv * (1 + S(n, w)))))
lauf.pruefe("ε_n' = (ε_n/x)(1 + S_n), ω_n' = ω_n/(x(1 + S_n(ω_n))) für n = 1…4", max(d) < mp.mpf(10) ** -20, mp.nstr(max(d), 3))

# ------------------------------------------------------------------ 5. Grenzidentität
x = mp.mpf("0.5")
ns = [5, 10, 20, 40, 80]
endlich = [S(n, omega_klein(n, x)) for n in ns]
y = x * mp.exp(-x)
Ns = [5, 10, 20, 40, 80]
video = [mp.fsum(mp.fprod(eps(m, y) for m in range(k - 1, N + 1)) for k in range(1, N + 1)) for N in Ns]
daten["grenz"] = {"n": ns, "endlich": [mp.nstr(v, 25) for v in endlich], "N": Ns, "video": [mp.nstr(v, 25) for v in video]}
lauf.pruefe("S_n(ω_n(1/2)) → x/(1−x) = 1 (n = 80: Abstand < 1e-20)", abs(endlich[-1] - 1) < mp.mpf(10) ** -20,
            {"n=5": mp.nstr(endlich[0], 12), "n=20": mp.nstr(endlich[2], 12), "n=80": mp.nstr(endlich[-1], 25)})
lauf.pruefe("Fassung aus dem Video (ω_∞ vorher eingesetzt) geht gegen x = 1/2, nicht gegen 1", abs(video[-1] - x) < mp.mpf(10) ** -15,
            {"N=80": mp.nstr(video[-1], 20)})

lauf.daten = daten
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogy(ns, [float(abs(v - 1)) for v in endlich], "o-", label="|S_n(ω_n(x)) − x/(1−x)|  (richtig: → 0)")
ax.semilogy(Ns, [float(abs(v - 1)) for v in video], "s--", label="|Σ Π ε_m(x e^{−x}) − x/(1−x)|  (Video: → 1/2)")
ax.set_xlabel("n bzw. N"); ax.set_title("Grenzidentität bei x = 1/2"); ax.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = ("Superwurzel-Formel exp(ω_{n−1}(ln Y)), h(x) = 1/ssrt(1/x), geschlossene Basisableitung und Ableitungen von ε_n, ω_n auf "
                        "≥ 20 Stellen bestätigt; die Grenzidentität gilt nur, wenn n → ∞ nach dem Einsetzen von ω_n geht.")
lauf.schreibe()
