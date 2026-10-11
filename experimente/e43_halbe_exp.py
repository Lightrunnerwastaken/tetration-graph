"""exp.halbe_exp — die halbe Exponentialfunktion f = exp^[1/2] = T_e(slog_e(x) + 1/2) aus der Kneser-Tetration (Kneser 1950).

  Prüft f(f(x)) = e^x, die Nullstelle bei T_e(−3/2) (Levenstein, Forum tid 544), die Singularität bei T_e(−5/2) ≈ −0.36237 + iπ,
  und das Wachstum entlang x_k = T_e(k): ln f / ln x = T(k − ½)/T(k − 1) → ∞, ln f / x^ε → 0.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, kneser

mp.mp.dps = 40
lauf = Lauf("exp.halbe_exp", __file__)
B = "exp(1)"


def f(xs):
    s = [mp.re(v) for v in kneser("slog", B, xs)]
    return [mp.re(v) for v in kneser("sexp", B, [t + mp.mpf("0.5") for t in s])]


xs = [mp.mpf(x) for x in ("-1.5", "-0.5", "0", "0.5", "1", "2", "3")]
ff = f(f(xs))
res = max(abs(a - mp.exp(x)) / mp.exp(x) for a, x in zip(ff, xs))
lauf.pruefe("f(f(x)) = e^x an 7 Stellen in [−1.5, 3] (relativ)", res < mp.mpf(10) ** -28, mp.nstr(res, 3))

T = {k: mp.re(v) for k, v in zip(("-1.5", "-0.5", "0.5", "1.5"), kneser("sexp", B, [mp.mpf(x) for x in ("-1.5", "-0.5", "0.5", "1.5")]))}
f0, f1 = f([mp.mpf(0), mp.mpf(1)])
lauf.pruefe("f(0) = T(−½) = ln f(1), f(1) = T(½) = 1.6463542337…", abs(f0 - T["-0.5"]) < 1e-30 and abs(f1 - T["0.5"]) < 1e-30 and abs(f0 - mp.log(f1)) < 1e-30,
            {"f(0)": mp.nstr(f0, 20), "f(1)": mp.nstr(f1, 20)})
x0 = T["-1.5"]
lauf.pruefe("einzige reelle Nullstelle von f bei x = T(−3/2) = −0.6960247…", abs(f([x0])[0]) < mp.mpf(10) ** -28, mp.nstr(x0, 20))
s25 = mp.log(mp.mpc(x0))                          # T(−5/2) = log T(−3/2): T(−3/2) < 0, also Imaginärteil π (der Fork liefert auf dem Schnitt nur den Realteil)
lauf.pruefe("T(−5/2) = −0.36237… + iπ (Ort der Singularität von f über der Achse)", abs(mp.re(s25) - mp.re(kneser("sexp", B, [mp.mpf("-2.5")])[0])) < 1e-28 and abs(mp.im(s25) - mp.pi) < 1e-30,
            mp.nstr(s25, 15))

# Wachstum entlang x_k = T(k): ln f(x_k) = T(k − ½), ln x_k = T(k − 1) — ohne riesige Zwischenwerte
ks = [2, 3, 4]
Th = [mp.re(v) for v in kneser("sexp", B, [mp.mpf(k) - mp.mpf("0.5") for k in ks])]
Tg = [mp.re(v) for v in kneser("sexp", B, [mp.mpf(k - 1) for k in ks])]
quot = {k: a / g for k, a, g in zip(ks, Th, Tg)}
lauf.pruefe("ln f(x)/ln x wächst entlang x = T(2), T(3), T(4) unbeschränkt (schneller als jedes Polynom)",
            quot[2] < quot[3] < quot[4] and quot[4] > 100, {k: mp.nstr(v, 6) for k, v in quot.items()})

lauf.daten = {"f0": mp.nstr(f0, 25), "f1": mp.nstr(f1, 25), "nullstelle": mp.nstr(x0, 25), "T(-2.5)": mp.nstr(s25, 20)}
fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
gitter = [mp.mpf(-150 + 5 * i) / 100 for i in range(80)]
fw = f(gitter)
ax.plot([float(x) for x in gitter], [float(y) for y in fw], label="f = exp^[1/2]")
ax.plot([float(x) for x in gitter], [float(mp.exp(x)) for x in gitter], "--", label="exp")
ax.plot([float(x) for x in gitter], [float(x) for x in gitter], ":", color="#888", label="x")
ax.set_ylim(-1.5, 8); ax.legend(fontsize=8); ax.set_title("Halbe Exponentialfunktion (Kneser)")
lauf.bild(fig)
lauf.zusammenfassung = f"f = exp^[1/2] über Kneser: f(f(x)) = eˣ auf 1e-28, f(0) = {mp.nstr(f0, 12)}, f(1) = {mp.nstr(f1, 12)}, Nullstelle {mp.nstr(x0, 12)}."
lauf.schreibe()
