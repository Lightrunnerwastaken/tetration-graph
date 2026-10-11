"""exp.slog — slog als geglätteter iterierter Logarithmus; Ableitungsidentitäten; Stellenzahl von Türmen.

  ⌈slog_a(y)⌉ = log*_a(y) für y > 1        (log* = Zahl der log_a bis ≤ 1)
  slog'(y) = ln a · a^y · slog'(a^y),  slog'(0) = ln a / T'(0)
  Stellen(a↑↑(n+1)) = ⌊a↑↑n · log10 a⌋ + 1
"""
from math import floor

import mpmath as mp

from gemeinsam import Lauf, cache_speichern, cauchy_ableitungen, kneser, plt

mp.mp.dps = 40
lauf = Lauf("exp.slog", __file__)
A = lambda ys: kneser("slog", "exp(1)", ys, dps=40)
T = lambda xs: kneser("sexp", "exp(1)", xs, dps=40)


def logstern(y):
    n = 0
    while y > 1:
        y = mp.log(y)
        n += 1
    return n


ys = [mp.mpf(v) for v in ("1.01", "1.5", "2.7", "3", "15", "16", "1e3", "1e6", "3814279", "1e100", "1e1000")]
werte = A(ys)
fehl = [(mp.nstr(y, 5), mp.nstr(mp.re(w), 8), logstern(y)) for y, w in zip(ys, werte) if mp.ceil(mp.re(w)) != logstern(y)]
lauf.pruefe("⌈slog_e(y)⌉ = log*(y) an 11 Punkten bis 10^1000", not fehl, fehl or "alle gleich")

d0 = cauchy_ableitungen(A, mp.mpf(0), mp.mpf("0.3"), 1, 64)[1]
Tstrich0 = cauchy_ableitungen(T, mp.mpf(0), mp.mpf("0.5"), 1, 64)[1]
lauf.pruefe("slog'(0) = ln e / T'(0)", abs(d0 - 1 / Tstrich0) < mp.mpf(10) ** -25, mp.nstr(d0, 20))
res = 0
for y in map(mp.mpf, ["-0.5", "0", "0.4"]):
    d_y = cauchy_ableitungen(A, y, mp.mpf("0.25"), 1, 64)[1]
    d_ey = cauchy_ableitungen(A, mp.e ** y, mp.mpf("0.25"), 1, 64)[1]
    res = max(res, abs(d_y - mp.e ** y * d_ey))
lauf.pruefe("slog'(y) = ln a · a^y · slog'(a^y) an 3 Stellen", res < mp.mpf(10) ** -25, mp.nstr(res, 3))

# Stellenzahl: Formel gegen exakte ganze Zahlen
exakt_ok = True
for a, n in [(2, 3), (2, 4), (3, 2), (3, 3), (4, 3), (5, 3), (7, 2)]:
    t = a
    for _ in range(n - 1):
        t = a ** t
    t_klein = a
    for _ in range(n - 2):
        t_klein = a ** t_klein
    formel = floor(t_klein * mp.log10(a)) + 1
    exakt_ok &= (len(str(t)) == formel)
lauf.pruefe("Stellenzahl ⌊a↑↑n·log10 a⌋+1 = exakt (7 Fälle bis 5↑↑3)", exakt_ok)
stellen_2_5 = floor(65536 * mp.log10(2)) + 1
mp.mp.dps = 60
stellen_3_4 = int(mp.floor(7625597484987 * mp.log10(3))) + 1
mp.mp.dps = 40
lauf.pruefe("2↑↑5 = 2^65536 hat 19729 Stellen", stellen_2_5 == 19729, stellen_2_5)
lauf.daten = {"slog_werte": {mp.nstr(y, 6): mp.nstr(mp.re(w), 15) for y, w in zip(ys, werte)}, "slog_strich_0": mp.nstr(d0, 25),
              "stellen_2_5": stellen_2_5, "stellen_3_4": stellen_3_4}

# Bild: slog gegen die Treppe log*
xs_log = [mp.mpf(10) ** (mp.mpf(i) / 40 - 1) for i in range(0, 241)]          # y von 0.1 bis 1e5
sl = A(xs_log)
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.semilogx([float(y) for y in xs_log], [float(mp.re(s)) for s in sl], lw=2, label="slog_e(y)")
ax.semilogx([float(y) for y in xs_log], [logstern(y) - (0 if y > 1 else 1) for y in xs_log], drawstyle="steps-post", label="log*(y) (Treppe)")
for k, v in [(0, 1), (1, mp.e), (2, mp.e ** mp.e), (3, mp.e ** mp.e ** mp.e)]:
    ax.axvline(float(v), ls=":", color="#aaa")
ax.set_xlabel("y (log. Skala)"); ax.set_ylabel("Höhe")
ax.set_title("slog ist die geglättete Treppe des iterierten Logarithmus")
ax.legend(fontsize=8)
lauf.bild(fig)
cache_speichern()
lauf.zusammenfassung = (f"⌈slog_e⌉ = log* an allen Punkten bis 10^1000; slog_e'(0) = {mp.nstr(d0, 16)} = 1/T_e'(0); "
                        f"Stellenformel exakt; 2↑↑5 hat 19729, 3↑↑4 hat {stellen_3_4} Stellen.")
lauf.schreibe()
