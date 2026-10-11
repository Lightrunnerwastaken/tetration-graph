"""exp.ln_basis — ln b als Ableitungsverhältnis und die Ableitung der Tetration nach der Basis.

  Aus ln T(x+1) = ln b · T(x):     (ln T)^(k)(x+1) = ln b · T^(k)(x),  also  ln b = T'(0)/T'(−1) = slog'(0)/slog'(1)
                                    und ln b = (T''(0) − T'(0)²)/T''(−1)   (b₇ aus TPID 7 ⇔ Wendepunkt bei −1)
  Basisableitung D = ∂_b T_b:       D(x+1) = ln b·T(x+1)·D(x) + T(x+1)·T(x)/b,  D(0) = 0
                                    D = T'·u  mit  u(x+1) − u(x) = T(x)/(b ln b T'(x))
  ganze Höhen (ohne Interpolation): ∂_b (b↑↑n) = Σ_{k<n} T(k)·(ln b)^(n−k−1)·∏_{j=k+1..n} T(j) / b
  slog-Seite G = ∂_b slog_b:        G(b^y) = G(y) − y·slog'(y)/(b ln b),  G(0) = G(1) = 0,  G(b) = −slog'(b)
 Geprüft an Kneser (b = 2, e, 10) über den Fork; Basisableitung per Richardson-Differenzen in b.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, kneser, cauchy_ableitungen
from stufen import basis_text

mp.mp.dps = 40
lauf = Lauf("exp.ln_basis", __file__)
BASEN = {"2": mp.mpf(2), "exp(1)": mp.e, "10": mp.mpf(10)}


def T_abl(basis, x0, r, m=3, N=128):
    return cauchy_ableitungen(lambda zs: kneser("sexp", basis, zs), mp.mpf(x0), mp.mpf(r), m_max=m, N=N)


def A_abl(basis, y0, r, m=2, N=128):
    return cauchy_ableitungen(lambda zs: kneser("slog", basis, zs), mp.mpf(y0), mp.mpf(r), m_max=m, N=N)


# 1. ln b als Verhältnis — an jeder Lösung gleich
tab, fehler1, fehler2, fehler3 = {}, 0, 0, 0
for name, b in BASEN.items():
    T0 = T_abl(name, 0, "0.5")
    Tm = T_abl(name, -1, "0.5")
    A0 = A_abl(name, 0, "0.4")
    A1 = A_abl(name, 1, "0.4")
    lb = mp.log(b)
    r1 = T0[1] / Tm[1]
    r2 = A0[1] / A1[1]
    r3 = (T0[2] - T0[1] ** 2) / Tm[2]
    r4 = A0[2] / (A1[1] + A1[2])
    fehler1 = max(fehler1, abs(r1 / lb - 1), abs(r2 / lb - 1))
    fehler2 = max(fehler2, abs(r3 / lb - 1), abs(r4 / lb ** 2 - 1))
    # dritte Ordnung: (ln T)'''(0) = ln b · T'''(−1)
    lnT3 = T0[3] - 3 * T0[2] * T0[1] + 2 * T0[1] ** 3
    fehler3 = max(fehler3, abs(lnT3 / (lb * Tm[3]) - 1))
    tab[name] = {"T'(0)": mp.nstr(T0[1], 15), "T'(-1)": mp.nstr(Tm[1], 15), "slog'(0)": mp.nstr(A0[1], 15),
                 "slog'(1)": mp.nstr(A1[1], 15), "T''(-1)": mp.nstr(Tm[2], 12)}
lauf.pruefe("ln b = T'(0)/T'(−1) = slog'(0)/slog'(1) (Kneser, b = 2, e, 10)", fehler1 < mp.mpf(10) ** -25, mp.nstr(fehler1, 3))
lauf.pruefe("ln b = (T''(0) − T'(0)²)/T''(−1) und (ln b)² = slog''(0)/(slog'(1) + slog''(1))", fehler2 < mp.mpf(10) ** -22, mp.nstr(fehler2, 3))
lauf.pruefe("(ln T)'''(0) = ln b · T'''(−1)", fehler3 < mp.mpf(10) ** -20, mp.nstr(fehler3, 3))

# Unabhängigkeit von der Lösung: T~(x) = T(x + ε sin 2πx) hat andere Steigungen, aber dasselbe Verhältnis
eps, h = mp.mpf("0.05"), mp.mpf(10) ** -12


def T_stoer(xs):
    return [mp.re(v) for v in kneser("sexp", "exp(1)", [x + eps * mp.sin(2 * mp.pi * x) for x in xs])]


v = T_stoer([h, -h, -1 + h, -1 - h])
s0, sm = (v[0] - v[1]) / (2 * h), (v[2] - v[3]) / (2 * h)
lauf.pruefe("gestörte Lösung T(x + 0.05 sin 2πx): T~'(0) ≠ T'(0), aber T~'(0)/T~'(−1) = 1 = ln e",
            abs(s0 / sm - 1) < mp.mpf(10) ** -15 and abs(s0 - mp.mpf(tab["exp(1)"]["T'(0)"])) > mp.mpf("0.1"),
            {"T~'(0)": mp.nstr(s0, 12), "Verhältnis": mp.nstr(s0 / sm, 18)})

# b₇ aus TPID 7: T''(−1) = 0
b7 = mp.mpf("10.22106341578140875546872955")
Tm7 = T_abl(basis_text(b7), -1, "0.5")
lauf.pruefe("bei b₇ = 10.2210634157814… (TPID 7) hat T einen Wendepunkt genau bei −1", abs(Tm7[2]) < mp.mpf(10) ** -20, mp.nstr(Tm7[2], 3))


# 2. Basisableitung bei ganzer Höhe: geschlossene Summe = Rekursion
def turm_werte(b, n):
    t = [mp.mpf(1)]
    for _ in range(n):
        t.append(mp.power(b, t[-1]))
    return t


def basis_summe(b, n):
    t, lb = turm_werte(b, n), mp.log(b)
    return mp.fsum(t[k] * lb ** (n - k - 1) * mp.fprod(t[k + 1:n + 1]) for k in range(n)) / b


def basis_rekursion(b, n):
    t, d = mp.mpf(1), mp.mpf(0)
    for _ in range(n):
        t_neu = mp.power(b, t)
        d = t_neu * (t / b + mp.log(b) * d)
        t = t_neu
    return d


res = max(abs(basis_summe(b, n) / basis_rekursion(b, n) - 1) for b in (mp.mpf("1.3"), mp.sqrt(2), mp.mpf(2), mp.e) for n in range(1, 5))
lauf.pruefe("∂_b(b↑↑n) = Σ T(k)(ln b)^(n−k−1)∏T(j)/b gleich der Rekursion (4 Basen, n = 1…4)", res < mp.mpf(10) ** -35, mp.nstr(res, 3))


# 3. Basisableitung bei gebrochener Höhe (Kneser, Richardson in b)
def D(b, xs, hb=mp.mpf(10) ** -7):
    def stufe(hh):
        p = kneser("sexp", basis_text(b + hh), xs)
        m = kneser("sexp", basis_text(b - hh), xs)
        return [mp.re(u - w) / (2 * hh) for u, w in zip(p, m)]
    d1, d2 = stufe(hb), stufe(hb / 2)
    return [(4 * y - x) / 3 for x, y in zip(d1, d2)]


xs = [mp.mpf("-0.5"), mp.mpf("0.25"), mp.mpf("0.5")]
res_d, res_u, werte_d = 0, 0, {}
for b in (mp.mpf(2), mp.e):
    bt = basis_text(b)
    pkt = xs + [x + 1 for x in xs] + [mp.mpf(0)]
    Dv = D(b, pkt)
    Tv = [mp.re(w) for w in kneser("sexp", bt, pkt)]
    lb = mp.log(b)
    n = len(xs)
    for i in range(n):
        rechts = lb * Tv[n + i] * Dv[i] + Tv[n + i] * Tv[i] / b
        res_d = max(res_d, abs(Dv[n + i] - rechts) / abs(Dv[n + i]))
        T1 = T_abl(bt, xs[i], "0.3", m=1, N=64)[1]
        T2 = T_abl(bt, xs[i] + 1, "0.3", m=1, N=64)[1]
        res_u = max(res_u, abs((Dv[n + i] / T2 - Dv[i] / T1) / (Tv[i] / (b * lb * T1)) - 1))
    werte_d[mp.nstr(b, 6)] = {mp.nstr(x, 3): mp.nstr(d, 12) for x, d in zip(pkt, Dv)}
    res_d = max(res_d, abs(Dv[-1]))                                   # D(0) = 0
lauf.pruefe("D(x+1) = ln b·T(x+1)·D(x) + T(x+1)T(x)/b und D(0) = 0 (Kneser b = 2, e; x = −½, ¼, ½)", res_d < mp.mpf(10) ** -12, mp.nstr(res_d, 3))
lauf.pruefe("u = D/T′ hat den Zuwachs T(x)/(b ln b T′(x))", res_u < mp.mpf(10) ** -11, mp.nstr(res_u, 3))

# 4. slog-Seite bei b = 2
b, hb = mp.mpf(2), mp.mpf(10) ** -7
ys = [mp.mpf(0), mp.mpf(1), mp.mpf("0.5"), mp.power(2, mp.mpf("0.5")), mp.mpf(2)]


def G(hh):
    p = kneser("slog", basis_text(b + hh), ys)
    m = kneser("slog", basis_text(b - hh), ys)
    return [mp.re(u - w) / (2 * hh) for u, w in zip(p, m)]


g1, g2 = G(hb), G(hb / 2)
Gv = [(4 * y - x) / 3 for x, y in zip(g1, g2)]
A05 = A_abl("2", mp.mpf("0.5"), "0.3", m=1, N=64)[1]
A2 = A_abl("2", mp.mpf(2), "0.3", m=1, N=64)[1]
r_slog = max(abs(Gv[0]), abs(Gv[1]), abs(Gv[3] - (Gv[2] - mp.mpf("0.5") * A05 / (2 * mp.log(2)))), abs(Gv[4] + A2))
lauf.pruefe("∂_b slog: G(0) = G(1) = 0, G(b^y) = G(y) − y slog'(y)/(b ln b), G(b) = −slog'(b) (b = 2)", r_slog < mp.mpf(10) ** -12, mp.nstr(r_slog, 3))

lauf.daten = {"ableitungen": tab, "basisableitung": werte_d, "G_b2": [mp.nstr(g, 12) for g in Gv],
              "T''(-1) bei b7": mp.nstr(Tm7[2], 5)}

# Bild: ln b aus Steigungen und D(x) über [−1, 2] für b = 2
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
gitter = [mp.mpf(-90 + 5 * i) / 100 for i in range(59)]
Dg = D(mp.mpf(2), gitter)
Tg = [mp.re(w) for w in kneser("sexp", "2", gitter)]
ax[0].plot([float(x) for x in gitter], [float(d) for d in Dg], label="∂_b T_b(x) bei b = 2")
ax[0].plot([float(x) for x in gitter], [float(t) for t in Tg], "--", label="T_2(x)")
ax[0].set_ylim(-1, 8); ax[0].set_xlabel("x"); ax[0].legend(fontsize=8)
ax[0].set_title("Ableitung nach der Basis (Kneser)")
namen = list(BASEN)
ax[1].bar(range(3), [float(mp.log(BASEN[n])) for n in namen], color="#ccd", label="ln b")
ax[1].plot(range(3), [float(mp.mpf(tab[n]["slog'(0)"]) / mp.mpf(tab[n]["slog'(1)"])) for n in namen], "o", label="slog′(0)/slog′(1)")
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["2", "e", "10"]); ax[1].legend(fontsize=8)
ax[1].set_title("ln b aus zwei Steigungen")
lauf.bild(fig)
lauf.zusammenfassung = ("ln b = T′(0)/T′(−1) = slog′(0)/slog′(1) und die Formen zweiter/dritter Ordnung an Kneser b = 2, e, 10 bestätigt, auch an einer gestörten Lösung; "
                        "b₇ ⇔ Wendepunkt bei −1; Basisableitung: Summenformel ganzer Höhen und die lineare Differenzengleichung bei gebrochener Höhe bestätigt.")
lauf.schreibe()
