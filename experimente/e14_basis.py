"""exp.basisableitung — Ableitung nach der Basis bei ganzer Höhe und am unendlichen Turm.

  t_1 = a, t_{n+1} = a^{t_n}:   t'_{n+1} = t_{n+1} (t_n / a + ln a · t'_n),  t'_1 = 1
  h(a) = a^h (Turm):            h'(a) = h² / (a (1 − h ln a))
  bei η:                         h(a) = e − √(2e³(η−a)/η) + O(η − a)
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 50
lauf = Lauf("exp.basisableitung", __file__)
eta = mp.e ** (1 / mp.e)


def turm(a, n):
    t = a
    for _ in range(n - 1):
        t = mp.power(a, t)
    return t


def turm_abl(a, n):
    t, d = a, mp.mpf(1)
    for _ in range(n - 1):
        t_neu = mp.power(a, t)
        d = t_neu * (t / a + mp.log(a) * d)
        t = t_neu
    return d


res = 0
for a in [mp.mpf("1.3"), mp.sqrt(2), mp.mpf("1.5"), mp.mpf(2), mp.e]:
    for n in (1, 2, 3, 4):
        if a > 1.6 and n == 4 and a != 2:
            continue
        num = mp.diff(lambda x: turm(x, n), a)
        res = max(res, abs(turm_abl(a, n) - num) / abs(num))
lauf.pruefe("Rekursion für d/da (a↑↑n) = numerische Ableitung (relativ)", res < mp.mpf(10) ** -35, mp.nstr(res, 3))
lauf.pruefe("Beispiel d/da a^a = a^a (1 + ln a) bei a = 2: 4(1 + ln 2)", abs(turm_abl(mp.mpf(2), 2) - 4 * (1 + mp.log(2))) < 1e-40,
            mp.nstr(turm_abl(mp.mpf(2), 2), 15))


def h(a):
    return -mp.lambertw(-mp.log(a)) / mp.log(a)


res = 0
for a in [mp.mpf("0.2"), mp.mpf("0.5"), mp.mpf("1.2"), mp.mpf("1.4"), mp.mpf("1.44")]:
    hh = mp.re(h(a))
    formel = hh ** 2 / (a * (1 - hh * mp.log(a)))
    res = max(res, abs(formel - mp.diff(lambda x: mp.re(h(x)), a)) / abs(formel))
lauf.pruefe("h'(a) = h²/(a(1 − h ln a)) am unendlichen Turm (5 Basen, relativ)", res < mp.mpf(10) ** -35, mp.nstr(res, 3))

abst = [mp.mpf(10) ** (-k) for k in range(2, 13)]
quot = [(mp.e - mp.re(h(eta - d))) / mp.sqrt(2 * mp.e ** 3 * d / eta) for d in abst]
lauf.pruefe("Wurzelgesetz: (e − h(a)) / √(2e³(η−a)/η) → 1", abs(quot[-1] - 1) < mp.mpf(10) ** -5, [mp.nstr(q, 10) for q in quot[-3:]])
korr = [(q - 1) / mp.sqrt(d) for q, d in zip(quot, abst)]
lauf.pruefe("Abweichung vom Wurzelgesetz fällt wie √(η − a) (nächste Ordnung)", abs(korr[-1] - korr[-2]) < 1e-3, mp.nstr(korr[-1], 8))

lauf.daten = {"wurzel_quotient": [mp.nstr(q, 12) for q in quot], "naechste_ordnung": mp.nstr(korr[-1], 10)}

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
aa = [mp.e ** (-mp.e) + (eta - mp.e ** (-mp.e)) * i / 300 for i in range(1, 300)]
ax[0].plot([float(a) for a in aa], [float(mp.re(h(a))) for a in aa], label="h(a) = a^(a^(a^…))")
nah = [eta - mp.mpf(10) ** (-mp.mpf(i) / 20) for i in range(20, 140)]
ax[0].plot([float(a) for a in nah], [float(mp.e - mp.sqrt(2 * mp.e ** 3 * (eta - a) / eta)) for a in nah], "--", label="e − √(2e³(η−a)/η)")
ax[0].axvline(float(eta), ls=":", color="#888")
ax[0].set_xlabel("Basis a"); ax[0].set_title("Unendlicher Turm und seine Wurzel bei η")
ax[0].legend(fontsize=8)
ax[1].loglog([float(d) for d in abst], [float(mp.e - mp.re(h(eta - d))) for d in abst], "o-", label="e − h(η − δ)")
ax[1].loglog([float(d) for d in abst], [float(mp.sqrt(2 * mp.e ** 3 * d / eta)) for d in abst], "--", label="Steigung 1/2")
ax[1].set_xlabel("δ = η − a"); ax[1].set_title("Quadratwurzel-Verzweigung")
ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Ableitungsrekursion und h'(a)-Formel auf ≤ 1e-35 bestätigt; (e − h)/√(2e³(η−a)/η) = {mp.nstr(quot[-1], 12)} bei η−a = 1e-12; "
                        f"nächste Ordnung ≈ {mp.nstr(korr[-1], 6)}·√(η−a).")
lauf.schreibe()
