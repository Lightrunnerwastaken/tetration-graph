"""exp.zwei_sinh — die 2sinh-Tetration (tommy1729, Forum tid 424/1144): T(x) = lim ln^[n](2sinh^[n+x+k](1)) mit T(0) = 1.

  2sinh(y) = e^y − e^−y hat den abstossenden Fixpunkt 0 mit Multiplikator 2; 2sinh^[t] über die Koenigs-Funktion (Grenzwerte mit asinh).
  Für grosse y ist 2sinh(y) ≈ e^y, darum ist ln^[n]∘2sinh^[n] schon für n = 3 auf 40 Stellen stabil. Verglichen mit Kneser:
  θ(x) = slog_Kneser(T_2sinh(x)) − x ist 1-periodisch, klein, aber nicht null.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, kneser

mp.mp.dps = 50
lauf = Lauf("exp.zwei_sinh", __file__)


def sigma(y, m=220):
    for _ in range(m):
        y = mp.asinh(y / 2)
    return y * mp.mpf(2) ** m


def psi(s, m=220):
    a = s / mp.mpf(2) ** m
    for _ in range(m):
        a = 2 * mp.sinh(a)
    return a


S1 = sigma(mp.mpf(1))


def F(t, n=None):
    """n = None: 2sinh anwenden, bis y > 60, dann noch einmal (Fehler danach ~ e^(−2e^60)); sonst genau n Schritte."""
    y = psi(mp.mpf(2) ** t * S1)
    m = 0
    while (n is None and y < 60) or (n is not None and m < n):
        y = 2 * mp.sinh(y)
        m += 1
    if n is None:
        y = 2 * mp.sinh(y)
        m += 1
    for _ in range(m):
        y = mp.log(y)
    return y


lauf.pruefe("ln^[n]∘2sinh^[n] stabil: n = 3, 4 und die Abbruchregel gleich auf 40 Stellen", abs(F(mp.mpf(0), 3) - F(mp.mpf(0), 4)) < mp.mpf(10) ** -40 and abs(F(mp.mpf(0)) - F(mp.mpf(0), 4)) < mp.mpf(10) ** -40, mp.nstr(F(mp.mpf(0)), 20))
k = mp.findroot(lambda t: F(t) - 1, mp.mpf("0.14"))
T = lambda x: F(x + k)
res = max(abs(T(x + 1) - mp.exp(T(x))) for x in (mp.mpf("-0.7"), mp.mpf("0.2"), mp.mpf("0.5")))
lauf.pruefe("T(x+1) = e^T(x) und T(0) = 1", res < mp.mpf(10) ** -40, mp.nstr(res, 3))
lauf.pruefe("T(−½) = 0.498743364531671 (Levenstein, tid 1144)", abs(T(mp.mpf("-0.5")) - mp.mpf("0.498743364531671")) < mp.mpf(10) ** -14, mp.nstr(T(mp.mpf("-0.5")), 20))

xs = [mp.mpf(i) / 24 for i in range(25)]
werte = [T(x) for x in xs]
sl = [mp.re(v) for v in kneser("slog", "exp(1)", werte)]
theta = [s - x for s, x in zip(sl, xs)]
spanne = max(theta) - min(theta)
lauf.pruefe("θ = slog_Kneser(T_2sinh(x)) − x ist periodisch (θ(0) = θ(1) = 0) und nicht konstant: Spanne > 10⁻⁵",
            abs(theta[0]) < 1e-30 and abs(theta[-1]) < 1e-30 and spanne > mp.mpf("1e-5"), {"Spanne": mp.nstr(spanne, 6), "max": mp.nstr(max(theta), 6), "min": mp.nstr(min(theta), 6)})
Th = kneser("sexp", "exp(1)", [mp.mpf("0.5")])[0]
lauf.daten = {"T_2sinh(0.5)": mp.nstr(T(mp.mpf("0.5")), 25), "Kneser(0.5)": mp.nstr(mp.re(Th), 25), "k": mp.nstr(k, 25), "theta_spanne": mp.nstr(spanne, 8)}

fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.8))
ax.plot([float(x) for x in xs], [float(t) * 1e4 for t in theta], "o-", ms=3)
ax.set_xlabel("x"); ax.set_ylabel("θ(x) · 10⁴")
ax.set_title("2sinh-Tetration gegen Kneser: slog_K(T_2sinh(x)) − x")
lauf.bild(fig)
lauf.zusammenfassung = (f"2sinh-Tetration: T(½) = {mp.nstr(T(mp.mpf('0.5')), 15)} gegen Kneser {mp.nstr(mp.re(Th), 15)}; "
                        f"die Höhenverschiebung θ wackelt mit Spanne {mp.nstr(spanne, 4)}.")
lauf.schreibe()
