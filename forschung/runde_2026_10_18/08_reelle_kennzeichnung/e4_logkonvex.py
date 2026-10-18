"""E4: Kneser — Vorzeichen von (log T')'' = (T'''T' − T''^2)/T'^2  (gleichwertig: 1/slog' konvex), x in (−2, 2]."""
import mpmath as mp, rk_lib as R
R.init(40)
x = mp.mpf('-1.98'); werte = []
while x <= 2.0001:
    c = R.kneser(x, 4, erlaube_exp=True)
    f2 = (6 * c[3] * c[1] - 4 * c[2] ** 2) / c[1] ** 2
    werte.append((x, f2))
    x += mp.mpf('0.02')
neg = [(mp.nstr(a, 4), mp.nstr(b, 5)) for a, b in werte if b <= 0]
print("Stellen mit (log T')'' <= 0:", neg[:5], "...", neg[-5:], " Anzahl", len(neg), "von", len(werte))
m = min(werte, key=lambda t: t[1]); print("Minimum", mp.nstr(m[0], 4), mp.nstr(m[1], 8))
for a, b in werte[::10]: print(mp.nstr(a, 4), mp.nstr(b, 8))
