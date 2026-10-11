"""exp.shell_thron — Fläche des Shell–Thron-Gebiets und seine äussersten Punkte (Forumsfragen 2022).

 Gebiet: b(t) = exp(t·e^{−t}), |t| ≤ 1 (fp.satz.shell_thron). b ist auf |t| < 1 schlicht (t e^{−t} sternförmig, Bild in |w| < e < π,
 dort ist exp injektiv), also Fläche = ∬ |b′(t)|² dA = π Σ n a_n² mit b(t) = Σ a_n t^n (Flächensatz).
 Gegenprobe: Randintegral ½ ∮ Im(conj(b) db) mit Trapezregel (spektral genau).
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 40
lauf = Lauf("exp.shell_thron", __file__)


def b(t):
    return mp.exp(t * mp.exp(-t))


def db(t):
    return b(t) * (1 - t) * mp.exp(-t)


a = mp.taylor(b, 0, 80)
flaeche = mp.pi * mp.fsum(n * a[n] ** 2 for n in range(1, 81))
NR = 4000
rand = mp.fsum(mp.im(mp.conj(b(mp.expjpi(2 * mp.mpf(j) / NR))) * db(mp.expjpi(2 * mp.mpf(j) / NR)) * 1j * mp.expjpi(2 * mp.mpf(j) / NR))
               for j in range(NR)) * (2 * mp.pi / NR) / 2
lauf.pruefe("Fläche aus dem Flächensatz π Σ n a_n² = Randintegral (auf 30 Stellen)", abs(flaeche - rand) < mp.mpf(10) ** -30,
            mp.nstr(flaeche, 30))
lauf.pruefe("Reihe abgeschnitten bei n = 80: Restglied n·a_n² < 1e-55", 80 * a[80] ** 2 < mp.mpf(10) ** -55, mp.nstr(a[80], 3))
lauf.pruefe("Schlichtheit: Re(t f′/f) = Re(1 − t) > 0 auf |t| < 1 für f = t e^{−t}, und |t e^{−t}| < e < π",
            all(mp.re(1 - 0.999 * mp.expjpi(2 * mp.mpf(j) / 360)) > 0 for j in range(360))
            and max(abs(0.999 * mp.expjpi(2 * mp.mpf(j) / 360) * mp.exp(-0.999 * mp.expjpi(2 * mp.mpf(j) / 360))) for j in range(360)) < mp.e)


# äusserste Punkte auf dem Rand |t| = 1
def rand_b(th):
    return b(mp.expj(th))


def extrem(f, start):
    th = mp.findroot(lambda s: mp.diff(f, s), start)
    return th, rand_b(th)


grob = [(mp.im(rand_b(2 * mp.pi * j / 720)), 2 * mp.pi * j / 720) for j in range(720)]
th_im, b_im = extrem(lambda s: mp.im(rand_b(s)), max(grob)[1])
grob_re = [(mp.re(rand_b(2 * mp.pi * j / 720)), 2 * mp.pi * j / 720) for j in range(720)]
th_li, b_li = extrem(lambda s: mp.re(rand_b(s)), min(grob_re)[1])
th_re, b_re = extrem(lambda s: mp.re(rand_b(s)), max(grob_re)[1])
lauf.pruefe("grösster Imaginärteil im Gebiet (auf dem Rand, Maximumprinzip)", mp.im(b_im) > 0 and abs(mp.diff(lambda s: mp.im(rand_b(s)), th_im)) < mp.mpf(10) ** -25,
            f"b = {mp.nstr(b_im, 20)} bei t = e^(i·{mp.nstr(th_im, 15)})")
lauf.pruefe("weder der rechteste noch der linkeste Punkt liegt auf der reellen Achse (dort nur [e^{−e}, η])",
            mp.re(b_re) > mp.e ** (1 / mp.e) + 0.5 and mp.re(b_li) < mp.exp(-mp.e) - 0.3 and abs(mp.im(b_re)) > 0.5 and abs(mp.im(b_li)) > 0.5,
            f"rechts: {mp.nstr(b_re, 15)}, links: {mp.nstr(b_li, 15)}")

# Vergleich: A. Robbins, Tetration Forum tid 263, 06.04.2009 (Randintegral für die obere Hälfte, verdoppelt)
ROB = {"flaeche": mp.mpf("8.05093328093950962342519537425"), "im": mp.mpc("0.8152027425068848021", "2.0166402199122300356"),
       "links": mp.mpc("-0.380979728647791984", "0.8997001955459000918"), "rechts": mp.mpc("2.0477905274644031305", "0.842045503530840715")}
lauf.pruefe("Robbins (2009): Fläche auf 30 Stellen, drei Extrempunkte auf 18 Stellen gleich",
            abs(flaeche - ROB["flaeche"]) < mp.mpf(10) ** -28 and abs(b_im - ROB["im"]) < mp.mpf(10) ** -18
            and abs(b_li - ROB["links"]) < mp.mpf(10) ** -17 and abs(mp.mpc(mp.re(b_re), abs(mp.im(b_re))) - ROB["rechts"]) < mp.mpf(10) ** -18,
            mp.nstr(flaeche - ROB["flaeche"], 3))

lauf.daten = {"flaeche": mp.nstr(flaeche, 35), "max_im": [mp.nstr(mp.re(b_im), 25), mp.nstr(mp.im(b_im), 25)], "theta_max_im": mp.nstr(th_im, 25),
              "links": [mp.nstr(mp.re(b_li), 25), mp.nstr(mp.im(b_li), 25)], "theta_links": mp.nstr(th_li, 25),
              "rechts": [mp.nstr(mp.re(b_re), 25), mp.nstr(mp.im(b_re), 25)], "theta_rechts": mp.nstr(th_re, 25),
              "koeffizienten": [mp.nstr(a[n], 15) for n in range(8)]}

fig, ax = plt.subplots(figsize=(6, 4.6))
rb = [rand_b(2 * mp.pi * j / 1000) for j in range(1001)]
ax.fill([float(mp.re(v)) for v in rb], [float(mp.im(v)) for v in rb], color="#3d6fb6", alpha=0.15)
ax.plot([float(mp.re(v)) for v in rb], [float(mp.im(v)) for v in rb], color="#3d6fb6", lw=1.2, label=f"Rand, Fläche {mp.nstr(flaeche, 8)}")
for v, name in ((b_im, "max Im"), (mp.conj(b_im), ""), (b_li, "links"), (b_re, "rechts"), (mp.e ** (1 / mp.e), "η"), (mp.exp(-mp.e), "e^(−e)")):
    ax.plot([float(mp.re(v))], [float(mp.im(v))], "o", color="#b5523b", ms=4)
    if name:
        ax.annotate(name, (float(mp.re(v)), float(mp.im(v))), (5, 5), textcoords="offset points", fontsize=8)
ax.set_aspect("equal"); ax.set_xlabel("Re b"); ax.set_ylabel("Im b"); ax.set_title("Shell–Thron-Gebiet: Fläche und Extrempunkte"); ax.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Fläche = {mp.nstr(flaeche, 20)}; grösster Imaginärteil bei b = {mp.nstr(b_im, 12)}; "
                        f"linkester Punkt {mp.nstr(b_li, 12)}, rechtester {mp.nstr(b_re, 12)} (beide abseits der reellen Achse).")
lauf.schreibe()
