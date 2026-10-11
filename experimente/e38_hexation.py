"""exp.hexation — reguläre Hexation unter η₅, Struktur jeder Hexation, unendliche Hexation, Fixpunkte der Pentation zur Basis 2.

 Hexation H_b(x+1) = P_b(H_b(x)), H_b(0) = 1, P_b = Pentation aus exp.pentation.
 Für η < b < η₅ hat P_b einen anziehenden reellen Fixpunkt p₁; reguläre Hexation H(x) = Ψ(μ^x σ(1)) dort (hexation.py).
 Struktur: H(−1) = 0, H(−2) = −1, H(−3) = −2; P_b nimmt auf ℝ nur Werte > L_T an, also H(x_c) = L_T für ein x_c ∈ (−3, −2)
 und H → −∞ bei s_b = x_c − 1 ∈ (−4, −3).
 Basis 2 (> η₅): P_2 hat keinen reellen Fixpunkt; komplexe Fixpunkte p (anziehend) und q (abstossend).
"""
import mpmath as mp

from gemeinsam import Lauf, kneser, plt
from hexation import HexationRegulaer, Pentation

mp.mp.dps = 40
lauf = Lauf("exp.hexation", __file__)

BASEN = ("1.5", "1.6", "1.7")
R = {b: HexationRegulaer(b) for b in BASEN}
GANZ = [mp.mpf(k) for k in (0, 1, -1, -2, -3)]
werte = {}
for b, r in R.items():
    v = r.H(GANZ + [mp.mpf(2), mp.mpf("0.5"), mp.mpf("1.5"), mp.mpf("-0.5"), mp.mpf(25)])
    werte[b] = v
    print(f"  b = {b}: p₁ = {mp.nstr(r.p1, 20)}, μ = {mp.nstr(r.mu, 12)}, L_T = {mp.nstr(r.Pe.L, 12)}, H(½) = {mp.nstr(v[6], 18)}")

lauf.pruefe("ganze Höhen: H(0) = 1, H(1) = b, H(−1) = 0, H(−2) = −1, H(−3) = −2 (drei Basen)",
            all(abs(v[0] - 1) + abs(v[1] - mp.mpf(b)) + abs(v[2]) + abs(v[3] + 1) + abs(v[4] + 2) < mp.mpf(10) ** -24 for b, v in werte.items()))
lauf.pruefe("H(2) = P(b) und P(H(½)) = H(1½)",
            all(abs(R[b].Pe.P([mp.mpf(b)])[0] - v[5]) + abs(R[b].Pe.P([v[6]])[0] - v[7]) < mp.mpf(10) ** -24 for b, v in werte.items()))
lauf.pruefe("beschränkt: H(x) → p₁ für x → ∞ (H(25) − p₁ wie μ^25)",
            all(0 < r.p1 - werte[b][9] < 2 * abs(r.sigma1) * r.mu ** 25 for b, r in R.items()),
            {b: mp.nstr(r.p1 - werte[b][9], 3) for b, r in R.items()})
lauf.pruefe("anziehend: 0 < μ = P′(p₁) < 1", all(0 < r.mu < 1 for r in R.values()), {b: mp.nstr(r.mu, 10) for b, r in R.items()})
gx = [mp.mpf(k) / 4 for k in range(-11, 25)]
kurven = {b: r.H(gx) for b, r in R.items()}
lauf.pruefe("streng wachsend auf (−2.75, 6]", all(all(k[i] < k[i + 1] for i in range(len(k) - 1)) for k in kurven.values()))

# Singularität: H(x_c) = L_T mit x_c ∈ (−3, −2), dann s_b = x_c − 1
sing = {}
for b, r in R.items():
    xa, xb = mp.mpf("-2.9"), mp.mpf("-2.1")
    fa, fb = r.H([xa])[0] - r.Pe.L, r.H([xb])[0] - r.Pe.L
    for _ in range(60):
        xc = xb - fb * (xb - xa) / (fb - fa)
        xa, fa, xb = xb, fb, xc
        fb = r.H([xb])[0] - r.Pe.L
        if abs(xb - xa) < mp.mpf(10) ** -22:
            break
    sing[b] = xb - 1
nahe = {b: R[b].H([s + mp.mpf("1e-6")])[0] for b, s in sing.items()}
lauf.pruefe("Singularität s_b ∈ (−4, −3): dort geht H gegen −∞ (H(s_b + 10⁻⁶) < −5)",
            all(-4 < s < -3 for s in sing.values()) and all(v < -5 for v in nahe.values()),
            {b: f"s = {mp.nstr(s, 15)}, H(s + 1e-6) = {mp.nstr(nahe[b], 6)}" for b, s in sing.items()})

# unendliche Hexation: b↑↑↑↑n = P^n(1) → p₁ (b < η₅), divergiert für b > η₅
folge = [mp.mpf(1)]
for _ in range(25):
    folge.append(mp.re(R["1.7"].Pe.P([folge[-1]])[0]))
lauf.pruefe("1.7↑↑↑↑n = P^n(1) wächst monoton gegen p₁ = 1.7↑↑↑↑∞",
            all(folge[i] < folge[i + 1] for i in range(len(folge) - 1)) and abs(folge[-1] - R["1.7"].p1) < mp.mpf(10) ** -8,
            mp.nstr(folge[-1], 15))
ueber = {}
for b, xmax in (("1.74", 40), ("2", 25)):                # 2↑↑↑x für x > 2.5 sprengt den Motor (2↑↑↑3 = 65536)
    Pe = Pentation(b)
    xs = [mp.mpf(k) / 10 for k in range(-30, xmax + 1)]
    ueber[b] = (Pe, min(mp.re(v) - x for v, x in zip(Pe.P(xs), xs)))
lauf.pruefe("über η₅ (1.74 auf [−3, 4], 2 auf [−3, 2.5]): P(x) > x — kein reeller Fixpunkt, die unendliche Hexation divergiert",
            all(m > 0 for _, m in ueber.values()), {b: mp.nstr(m, 6) for b, (_, m) in ueber.items()})

# Basis 2: komplexe Fixpunkte der Pentation
Pe2 = ueber["2"][0]
p, ap = Pe2.fixpunkt(mp.mpc("1.7937303432391", "1.5736302435118"))
q, aq = Pe2.fixpunkt(mp.mpc("-1.9929484359323", "1.4451284966392"), r=mp.mpf("0.1"))
lauf.pruefe("Basis 2: Fixpunkt p = 1.79373 + 1.57363i ist anziehend (|μ| = 0.6087)", abs(ap[0] - p) < mp.mpf(10) ** -25 and abs(ap[1]) < 1,
            f"p = {mp.nstr(p, 20)}, μ = {mp.nstr(ap[1], 12)}")
lauf.pruefe("Basis 2: Fixpunkt q = −1.99295 + 1.44513i ist abstossend (|μ| = 2.84)", abs(aq[0] - q) < mp.mpf(10) ** -25 and abs(aq[1]) > 1,
            f"q = {mp.nstr(q, 20)}, μ = {mp.nstr(aq[1], 12)}")

lauf.daten = {b: {"p1": mp.nstr(r.p1, 25), "mu": mp.nstr(r.mu, 20), "L_T": mp.nstr(r.Pe.L, 20), "sigma1": mp.nstr(r.sigma1, 20),
                  "H_halb": mp.nstr(werte[b][6], 25), "H_1halb": mp.nstr(werte[b][7], 25), "H_minus_halb": mp.nstr(werte[b][8], 25),
                  "H2": mp.nstr(werte[b][5], 25), "s_b": mp.nstr(sing[b], 20)} for b, r in R.items()}
lauf.daten["basis2"] = {"p": [mp.nstr(mp.re(p), 25), mp.nstr(mp.im(p), 25)], "mu_p": [mp.nstr(mp.re(ap[1]), 18), mp.nstr(mp.im(ap[1]), 18)],
                        "q": [mp.nstr(mp.re(q), 25), mp.nstr(mp.im(q), 25)], "mu_q": [mp.nstr(mp.re(aq[1]), 18), mp.nstr(mp.im(aq[1]), 18)]}
lauf.daten["folge_1.7"] = [mp.nstr(v, 15) for v in folge[:10]]

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for (b, k), f in zip(kurven.items(), ("#2f8f83", "#3d6fb6", "#b5523b")):
    ax[0].plot([float(x) for x in gx], [float(y) for y in k], color=f, label=f"b = {b}")
    ax[0].axhline(float(R[b].p1), ls=":", color=f, lw=0.8)
    ax[0].axvline(float(sing[b]), ls="--", color=f, lw=0.6)
ax[0].plot([-3, 3], [-3, 3], ":", color="#8a93a3", lw=0.8)
ax[0].set_ylim(-4, 3); ax[0].set_xlim(-4, 6); ax[0].set_xlabel("x"); ax[0].set_title("Reguläre Hexation b↑↑↑↑x (b < η₅)"); ax[0].legend(fontsize=8)
r7 = R["1.7"]
xx = [mp.mpf(k) / 10 for k in range(-28, 31)]
tx = [x for x in xx if x > mp.mpf("-1.95")]
ax[1].plot([float(x) for x in tx], [float(mp.re(v)) for v in kneser("sexp", "1.7", tx)], color="#8a93a3", label="Tetration 1.7↑↑x")
ax[1].plot([float(x) for x in xx], [float(v) for v in r7.Pe.P(xx)], color="#7a49c2", label="Pentation 1.7↑↑↑x")
hx = [x for x in xx if x > sing["1.7"] + mp.mpf("0.02")]
ax[1].plot([float(x) for x in hx], [float(v) for v in r7.H(hx)], color="#b5523b", label="Hexation 1.7↑↑↑↑x")
ax[1].plot([-1, 0, 1], [0, 1, 1.7], "o", color="#222", ms=3)
ax[1].set_ylim(-3.2, 4); ax[1].set_xlabel("x"); ax[1].set_title("Basis 1.7: drei Stufen durch (−1, 0), (0, 1), (1, b)"); ax[1].legend(fontsize=7.5)
lauf.bild(fig)
lauf.zusammenfassung = ("Reguläre Hexation unter η₅: " + ", ".join(f"{b}↑↑↑↑½ = {mp.nstr(werte[b][6], 12)}, {b}↑↑↑↑∞ = {mp.nstr(R[b].p1, 12)}" for b in BASEN)
                        + f"; H(−3) = −2, Singularität bei s ∈ (−4, −3); über η₅ kein reeller Fixpunkt der Pentation; Basis 2: p anziehend, q abstossend.")
lauf.schreibe()
