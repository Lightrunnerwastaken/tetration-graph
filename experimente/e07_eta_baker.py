"""exp.eta_baker — exp_η am Fixpunkt e ist linear konjugiert zu e^u − 1; formale halbe Iterierte divergiert.

(1) η^{e+eu} − e − e(e^u − 1) = 0 an Zufallspunkten (auch komplex).
(2) Koeffizienten c_n der formalen halben Iterierten g von f(u) = e^u − 1 (g∘g = f), exakt rational.
    Wachstum |c_n|: faktoriell (Gevrey-1) ⇒ Konvergenzradius 0, verträglich mit Bakers Satz.
"""
import random
from fractions import Fraction
from math import factorial, lgamma, log

import mpmath as mp

from gemeinsam import Lauf, plt

lauf = Lauf("exp.eta_baker", __file__)
mp.mp.dps = 60
eta = mp.e ** (1 / mp.e)
random.seed(2)
res = 0
for _ in range(20):
    u = mp.mpc(random.uniform(-2, 2), random.uniform(-2, 2))
    res = max(res, abs(mp.power(eta, mp.e + mp.e * u) - mp.e - mp.e * (mp.exp(u) - 1)))
lauf.pruefe("η^{e+eu} = e + e(e^u − 1) an 20 komplexen Zufallspunkten", res < mp.mpf(10) ** -50, mp.nstr(res, 3))

N = 60
f = [Fraction(0), Fraction(1)] + [Fraction(1, factorial(n)) for n in range(2, N + 1)]


def compose(a, b, n):
    """Koeffizienten von a(b(u)) bis Grad n (b[0] = 0)."""
    out = [Fraction(0)] * (n + 1)
    p = [Fraction(0)] * (n + 1)
    p[0] = Fraction(1)
    for k in range(1, n + 1):
        q = [Fraction(0)] * (n + 1)
        for i, pi in enumerate(p):
            if pi:
                for j in range(1, n + 1 - i):
                    if b[j]:
                        q[i + j] += pi * b[j]
        p = q
        if a[k]:
            for i in range(n + 1):
                out[i] += a[k] * p[i]
    return out


c = [Fraction(0), Fraction(1)] + [Fraction(0)] * (N - 1)
for n in range(2, N + 1):
    gg = compose(c, c, n)
    # [u^n] g(g(u)) = 2 c_n + (Rest ohne c_n); also c_n = (f_n − Rest)/2
    c[n] = (f[n] - gg[n]) / 2
gg = compose(c, c, N)
lauf.pruefe("g∘g = e^u − 1 bis Grad 60 exakt", all(gg[n] == f[n] for n in range(1, N + 1)))
lauf.pruefe("bekannte Anfangskoeffizienten c_2 = 1/4, c_3 = 1/48", c[2] == Fraction(1, 4) and c[3] == Fraction(1, 48),
            [str(c[2]), str(c[3])])
betrag = [abs(float(x)) for x in c]
# Wachstum: ln|c_n| − ln Γ(n) sollte linear in n verlaufen (Gevrey-1), Wurzeltest |c_n|^{1/n} wächst unbeschränkt
nz = [n for n in range(10, N + 1) if betrag[n] > 0]
wurzel = [betrag[n] ** (1 / n) for n in nz]
lauf.pruefe("|c_n|^{1/n} wächst (n = 30 → 60) — kein endlicher Konvergenzradius sichtbar",
            max(wurzel[-6:]) > 2 * max(wurzel[:6]), [round(w, 4) for w in (wurzel[0], wurzel[-1])])
gevrey = [log(betrag[n]) - lgamma(n) for n in nz]
lauf.daten = {"c_n": [str(x) for x in c[:16]], "wurzeltest": dict(zip(nz, wurzel)), "ln_c_minus_lnGamma": dict(zip(nz, gevrey))}

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].semilogy(nz, [betrag[n] for n in nz], "o", ms=3)
ax[0].set_xlabel("n")
ax[0].set_ylabel("|c_n|")
ax[0].set_title("Formale halbe Iterierte von e^u − 1")
ax[1].plot(nz, gevrey, "o", ms=3)
ax[1].set_xlabel("n")
ax[1].set_ylabel("ln|c_n| − ln Γ(n)")
ax[1].set_title("Gevrey-1-Test (lineare Hülle ⇒ faktorielles Wachstum)")
lauf.bild(fig)
lauf.zusammenfassung = (f"Konjugation exakt (Residuum {mp.nstr(res, 2)}). Halbe Iterierte: c_2 = 1/4, c_3 = 1/48, …; "
                        f"|c_n|^(1/n) wächst von {wurzel[0]:.3f} (n=10) auf {wurzel[-1]:.3f} (n=60) — divergent, wie Baker.")
lauf.schreibe()
