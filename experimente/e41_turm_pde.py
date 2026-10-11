"""exp.turm_pde — Kobis partielle Differentialgleichung für Türme (Forum tid 417, 2010) und die Basisableitung der Kneser-Tetration.

  f_n = exp_a^[n] (ganze n ≥ −1),  S_a(x) = (1/(a ln²a)) Σ_{k≥−1} 1/f_k'(x)   (Vorwärtsbahn; f_{−1} = log_a)
  Satz:  y = f_n(a, x) erfüllt  ∂_a y = S_a(x)·∂_x y − S_a(y)  für jedes ganze n (Identität).
  Für gebrochene Höhen y = α^{-1}(α(x) + h) gilt sie genau dann, wenn p := ∂_a α − S·α' konstant ist;
  p ist immer 1-periodisch in α. Gemessen: p für Kneser (b = 2, e) über eine Periode — fast konstant, mit kleinem Wackeln.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, kneser, cauchy_ableitungen
from stufen import basis_text

mp.mp.dps = 40
lauf = Lauf("exp.turm_pde", __file__)


def iteriert(a, x, n):
    """f_n(a, x) und ∂_x f_n für ganze n ≥ −1."""
    la = mp.log(a)
    if n == -1:
        return mp.log(x) / la, 1 / (x * la)
    f, d = x, mp.mpf(1)
    for _ in range(n):
        f = mp.power(a, f)
        d = la * f * d
    return f, d


def S(a, x, kmax=40):
    la = mp.log(a)
    tot = x * la + 1                           # k = −1 und k = 0
    f, d = x, mp.mpf(1)
    for _ in range(kmax):
        f = mp.power(a, f)
        d = la * f * d
        if abs(1 / d) < mp.mpf(10) ** -48:
            break
        tot += 1 / d
    return tot / (a * la ** 2)


# 1. Identität für ganze Höhen
res = 0
for a in (mp.mpf("1.7"), mp.mpf(2), mp.e, mp.mpf(3)):
    for x in (mp.mpf("0.3"), mp.mpf(1), mp.mpf("1.5")):
        for n in (1, 2, 3):
            y, dy = iteriert(a, x, n)
            da = mp.diff(lambda aa: iteriert(aa, x, n)[0], a)
            res = max(res, abs(da - (S(a, x) * dy - S(a, y))) / abs(da))
lauf.pruefe("∂_a f_n = S(x)·∂_x f_n − S(f_n) für n = 1, 2, 3 (4 Basen, 3 Stellen, relativ)", res < mp.mpf(10) ** -30, mp.nstr(res, 3))

# S erfüllt dieselbe inhomogene Gleichung wie γ = ∂_a slog / slog':  S(a^x) = a^x ln a·S(x) − x a^(x−1)
res = max(abs(S(mp.mpf(2), mp.power(2, x)) - (mp.power(2, x) * mp.log(2) * S(mp.mpf(2), x) - x * mp.power(2, x - 1))) for x in (mp.mpf("0.2"), mp.mpf("0.7"), mp.mpf(1)))
lauf.pruefe("S(a^x) = a^x ln a · S(x) − x·a^(x−1)", res < mp.mpf(10) ** -35, mp.nstr(res, 3))

# 2. p = ∂_b slog − S·slog' für Kneser über eine Periode y ∈ [0, 1] (slog von −1 bis 0)
M = 24
ergebnis, kurven = {}, {}
for name, b in (("2", mp.mpf(2)), ("exp(1)", mp.e)):
    ys = [mp.mpf(k) / M for k in range(M + 1)]
    h = mp.mpf(10) ** -7

    def G(hh):
        p_ = kneser("slog", basis_text(b + hh), ys)
        m_ = kneser("slog", basis_text(b - hh), ys)
        return [mp.re(u - w) / (2 * hh) for u, w in zip(p_, m_)]
    g1, g2 = G(h), G(h / 2)
    Gs = [(4 * y - x) / 3 for x, y in zip(g1, g2)]
    A = [mp.re(v) for v in kneser("slog", name, ys)]
    Ad = [cauchy_ableitungen(lambda zs: kneser("slog", name, zs), y, mp.mpf("0.25"), m_max=1, N=96)[1] for y in ys]
    p = [g - S(b, y) * ad for g, y, ad in zip(Gs, ys, Ad)]
    # Fourier in A: p als Funktion von A ∈ [−1, 0]; Mittel und erste Harmonische per Trapez auf ungleichem Gitter
    mittel = mp.fsum((p[i] + p[i + 1]) / 2 * (A[i + 1] - A[i]) for i in range(M))
    c1 = mp.fsum(((p[i] - mittel) * mp.expjpi(-2 * A[i]) + (p[i + 1] - mittel) * mp.expjpi(-2 * A[i + 1])) / 2 * (A[i + 1] - A[i]) for i in range(M))
    ergebnis[name] = {"mittel": mp.nstr(mittel, 12), "harm1": mp.nstr(abs(c1), 4), "spanne": mp.nstr(max(p) - min(p), 4),
                      "periodisch": mp.nstr(abs(p[0] - p[-1]), 3)}
    kurven[name] = (A, p, mittel)
lauf.pruefe("p ist 1-periodisch in slog: p(slog = −1) = p(slog = 0) (b = 2, e)", all(mp.mpf(e["periodisch"]) < mp.mpf(10) ** -15 for e in ergebnis.values()),
            {k: v["periodisch"] for k, v in ergebnis.items()})
lauf.pruefe("Kneser: p fast konstant, aber nicht konstant — Spanne klein gegen das Mittel, doch weit über der Rechengenauigkeit",
            all(mp.mpf("1e-8") < mp.mpf(e["spanne"]) < mp.mpf("1e-2") * abs(mp.mpf(e["mittel"])) for e in ergebnis.values()),
            ergebnis)

lauf.daten = {"p_kneser": ergebnis}
fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.8))
for name, (A, p, m) in kurven.items():
    ax.plot([float(a) for a in A], [float((q - m) * 10 ** 4) for q in p], "o-", ms=3, label=f"b = {'e' if name == 'exp(1)' else name}: Mittel {mp.nstr(m, 8)}")
ax.set_xlabel("slog_b(y)"); ax.set_ylabel("(p − Mittel) · 10⁴")
ax.set_title("Kneser weicht von Kobis Gleichung um ein kleines Wackeln ab")
ax.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Kobis Gleichung gilt für ganze Höhen exakt; für Kneser ist p = ∂_b slog − S·slog′ periodisch mit Mittel {ergebnis['2']['mittel']} (b = 2) "
                        f"und Spanne {ergebnis['2']['spanne']} — Kneser erfüllt sie bei gebrochener Höhe also nur fast.")
lauf.schreibe()
