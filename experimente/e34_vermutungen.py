"""exp.vermutungen — Stichproben zu offenen Vermutungen aus Forum und Wikis.

 TPID 7 (andydude 2009):  0 < T(y)·T(−y) < 1 für 0 < |y| < 1 (regulär unter η, Kneser darüber).
 Halbebene (Nixon 2017, Forum tid 1168): bildet die reguläre Tetration (1 < b < η) Re z ≥ 0 in Re w > 0 ab?
   Hier: aus c_k ≤ 0 (kl.verm.vollmonoton) folgt |T(z) − L| ≤ L − 1 auf Re z ≥ 0, also Re T ≥ 1. Gitterprobe.
 ∫_{−1}^0 T_b = 1/2: gilt für die lineare Näherung immer; für welche Basis bei Kneser?
 TPID 1 (bo198214 2008): Eigenwerte der abgeschnittenen Carleman-Matrix von b^x gehen gegen (ln a)^k = λ^k.
"""
import mpmath as mp

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, cauchy_ableitungen, gauss_legendre, kneser, plt

mp.mp.dps = 40
lauf = Lauf("exp.vermutungen", __file__)
YS = [mp.mpf(k) / 20 for k in range(1, 20)]


def regulaer(b):
    u = mp.log(b)
    L = mp.re(-mp.lambertw(-u) / u)
    F = AmFixpunkt(u, L, N=90)
    s1 = F.sigma(mp.mpf(1))
    return F, s1, L


# ------------------------------------------------------------------ TPID 7
produkte = {}
for b in (mp.mpf("1.1"), mp.mpf("1.3"), mp.sqrt(2), mp.mpf("1.44")):
    F, s1, L = regulaer(b)
    T = lambda z: mp.re(F.Psi(F.lam ** z * s1))
    produkte[mp.nstr(b, 4)] = [T(y) * T(-y) for y in YS]
for b in ("1.5", "2", "exp(1)", "3", "5", "10", "20", "100"):
    w = [mp.re(v) for v in kneser("sexp", b, YS + [-y for y in YS])]
    produkte[b] = [w[i] * w[i + len(YS)] for i in range(len(YS))]
gut = [k for k in produkte if k not in ("20", "100")]
lauf.pruefe("TPID 7: 0 < T(y)T(−y) < 1 und monoton fallend für y = 0.05 … 0.95 (vier reguläre Basen; Kneser 1.5 … 10)",
            all(0 < p < 1 for k in gut for p in produkte[k]) and all(produkte[k][i] > produkte[k][i + 1] for k in gut for i in range(len(YS) - 1)),
            {k: mp.nstr(produkte[k][9], 6) for k in gut})
lauf.pruefe("Gegenbeispiel zur Kneser-Fassung: Basis 20 und 100 haben T(½)T(−½) > 1",
            produkte["20"][9] > 1 and produkte["100"][9] > 1, {k: mp.nstr(produkte[k][9], 6) for k in ("20", "100")})

# lokal: T(y)T(−y) = exp((ln T)''(0) y² + …), (ln T)''(0) = T''(0) − T'(0)²; die Grenze b₇ ist die Nullstelle
def logkonkav(b):
    abl = cauchy_ableitungen(lambda zs: kneser("sexp", b if isinstance(b, str) else mp.nstr(b, 40), zs), mp.mpf(0), mp.mpf("0.3"), m_max=2, N=128)
    return abl[2] - abl[1] ** 2   # Radius 0.3: bei grossen Basen gab r = 1, N = 48 Abtastfehler ~1e-9 (b₇ falsch ab der 8. Stelle)


lauf.pruefe("lokal bei y → 0: 1 − T(y)T(−y) ≈ (T′(0)² − T″(0)) y² (Basis e: 0.6490)",
            abs((1 - produkte["exp(1)"][0]) / YS[0] ** 2 + logkonkav("exp(1)")) < mp.mpf("0.01"),
            f"{mp.nstr((1 - produkte['exp(1)'][0]) / YS[0] ** 2, 6)} gegen {mp.nstr(-logkonkav('exp(1)'), 6)}")
d0, d1 = mp.mpf(10), mp.mpf(12)
k0, k1 = logkonkav(d0), logkonkav(d1)
for _ in range(30):
    d2 = d1 - k1 * (d1 - d0) / (k1 - k0)
    d0, k0, d1 = d1, k1, d2
    k1 = logkonkav(d1)
    if abs(d1 - d0) < mp.mpf(10) ** -22:
        break
b7 = d1
probe = {}
for b in (b7 - mp.mpf("0.5"), b7 + mp.mpf("0.5")):
    w = [mp.re(v) for v in kneser("sexp", mp.nstr(b, 40), YS[:4] + [-y for y in YS[:4]])]
    probe[mp.nstr(b, 6)] = [w[i] * w[i + 4] for i in range(4)]
lauf.pruefe("Grenzbasis b₇ mit T″(0) = T′(0)² (Kneser); knapp darunter Produkt < 1, knapp darüber > 1 (y = 0.05 … 0.2)",
            10 < b7 < 20 and all(p < 1 for p in list(probe.values())[0]) and all(p > 1 for p in list(probe.values())[1]),
            f"b₇ = {mp.nstr(b7, 20)}")

# ------------------------------------------------------------------ Halbebene
halb = {}
for b in (mp.mpf("1.1"), mp.mpf("1.3"), mp.sqrt(2), mp.mpf("1.44")):
    F, s1, L = regulaer(b)
    per = 2 * mp.pi / abs(mp.log(F.lam))
    rmin, dmax = mp.inf, mp.mpf(0)
    for x in [mp.mpf(i) / 4 for i in range(0, 13)]:
        for j in range(-24, 25):
            v = F.Psi(F.lam ** mp.mpc(x, per * j / 48) * s1)
            rmin = min(rmin, mp.re(v))
            dmax = max(dmax, abs(v - L) / (L - 1))
    halb[mp.nstr(b, 4)] = (rmin, dmax)
lauf.pruefe("Halbebene: Re T(z) ≥ 1 und |T(z) − L| ≤ L − 1 auf Re z ∈ [0, 3] × eine volle Periode (je 637 Punkte)",
            all(r > 1 - mp.mpf(10) ** -25 and d < 1 + mp.mpf(10) ** -25 for r, d in halb.values()),
            {k: f"min Re T = {mp.nstr(r, 8)}, max |T−L|/(L−1) = {mp.nstr(d, 8)}" for k, (r, d) in halb.items()})

# ------------------------------------------------------------------ ∫_{−1}^0 T_b = 1/2
def integral(b):
    return gauss_legendre(lambda xs: [mp.re(v) for v in kneser("sexp", b if isinstance(b, str) else mp.nstr(b, 40), xs)], -1, 0, n=40)


I_e = integral("exp(1)")
lauf.pruefe("Kontrolle: ∫_{−1}^0 T_e = 0.49921818282999133972… (gr.beob.integrale_kneser)", abs(I_e - mp.mpf("0.49921818282999133972")) < mp.mpf(10) ** -18)
c0, c1 = mp.mpf(2), mp.mpf("2.5")
h0, h1 = integral(c0) - mp.mpf(1) / 2, integral(c1) - mp.mpf(1) / 2
werte_I = {"2": h0 + mp.mpf(1) / 2, "2.5": h1 + mp.mpf(1) / 2}
for _ in range(20):
    c2 = c1 - h1 * (c1 - c0) / (h1 - h0)
    c0, h0, c1 = c1, h1, c2
    h1 = integral(c1) - mp.mpf(1) / 2
    if abs(c1 - c0) < mp.mpf(10) ** -25:
        break
b_halb = c1
for b in ("1.5", "3", "10"):
    werte_I[b] = integral(b)
lauf.pruefe("∫_{−1}^0 T_b = 1/2 genau bei einer Basis zwischen 2 und e (Kneser)", 2 < b_halb < mp.e and abs(h1) < mp.mpf(10) ** -25,
            f"b = {mp.nstr(b_halb, 25)}")
lauf.pruefe("∫_{−1}^0 T_b fällt mit b (1.5, 2, 2.5, e, 3, 10)",
            werte_I["1.5"] > werte_I["2"] > werte_I["2.5"] > I_e > werte_I["3"] > werte_I["10"], {k: mp.nstr(v, 10) for k, v in werte_I.items()})

# ------------------------------------------------------------------ TPID 1: Carleman-Eigenwerte
def carleman_ew(b, N):
    lb = mp.log(b)
    M = mp.matrix(N, N)
    for j in range(N):
        for k in range(N):
            M[j, k] = (j * lb) ** k / mp.factorial(k)
    return sorted([mp.re(e) for e in mp.eig(M, left=False, right=False)], key=lambda z: -abs(z))


carl = {}
for b in (mp.mpf("1.3"), mp.sqrt(2)):
    lam = mp.re(-mp.lambertw(-mp.log(b)))
    carl[mp.nstr(b, 4)] = {N: max(abs(e - lam ** k) for k, e in enumerate(carleman_ew(b, N)[:4])) for N in (8, 16, 24, 32)}
lauf.pruefe("TPID 1: die ersten vier Eigenwerte der N×N-Carleman-Matrix nähern sich λ^k (Abstand fällt mit N)",
            all(v[8] > v[16] > v[24] > v[32] for v in carl.values()) and carl["1.3"][32] < mp.mpf(10) ** -9,
            {b: {N: mp.nstr(d, 3) for N, d in v.items()} for b, v in carl.items()})

lauf.daten = {"tpid7_y_halb": {k: mp.nstr(v[9], 15) for k, v in produkte.items()}, "b7": mp.nstr(b7, 25),
              "logkonkav_e": mp.nstr(logkonkav("exp(1)"), 20),
              "halbebene": {k: [mp.nstr(r, 15), mp.nstr(d, 15)] for k, (r, d) in halb.items()},
              "b_integral_halb": mp.nstr(b_halb, 30), "integrale": {k: mp.nstr(v, 20) for k, v in werte_I.items()},
              "carleman_abstand": {b: {N: mp.nstr(d, 4) for N, d in v.items()} for b, v in carl.items()}}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
farben = ["#9fb7d9", "#6f93c6", "#3d6fb6", "#244f8f", "#e8b8a6", "#e0a58f", "#cf7d61", "#b5523b", "#8f3a26", "#6b2a1b", "#7a49c2", "#4a1d12"]
for (k, v), f in zip(produkte.items(), farben):
    ax[0].plot([0] + [float(y) for y in YS] + [1], [1] + [float(p) for p in v] + [0], color=f, lw=1.2,
               ls="--" if k in ("20", "100") else "-", label=("b = " + k.replace("exp(1)", "e")))
ax[0].axhline(1, color="#8a93a3", lw=0.6)
ax[0].set_ylim(0, 1.7); ax[0].set_xlabel("y"); ax[0].set_ylabel("T(y)·T(−y)"); ax[0].set_title(f"TPID 7: unter 1 bis b₇ = {mp.nstr(b7, 5)}")
ax[0].legend(fontsize=6, ncol=2)
F, s1, L = regulaer(mp.sqrt(2))
th = [2 * mp.pi * j / 200 for j in range(201)]
ax[1].plot([float(L + (L - 1) * mp.cos(t)) for t in th], [float((L - 1) * mp.sin(t)) for t in th], ":", color="#8a93a3", label="|w − L| = L − 1")
per = 2 * mp.pi / abs(mp.log(F.lam))
for x, f in ((0, "#b5523b"), (mp.mpf("0.5"), "#7a49c2"), (mp.mpf(1), "#2f8f83"), (mp.mpf(2), "#3d6fb6")):
    vals = [F.Psi(F.lam ** mp.mpc(x, per * j / 200) * s1) for j in range(-100, 101)]
    ax[1].plot([float(mp.re(v)) for v in vals], [float(mp.im(v)) for v in vals], color=f, lw=1.1, label=f"T(x + iy), x = {mp.nstr(x, 2)}")
ax[1].axvline(1, color="#8a93a3", lw=0.6)
ax[1].set_aspect("equal"); ax[1].set_title("Basis √2: Bild der rechten Halbebene"); ax[1].legend(fontsize=7)
lauf.bild(fig)
lauf.zusammenfassung = (f"TPID 7 hält regulär und für Kneser bis b₇ = {mp.nstr(b7, 12)} (T″(0) = T′(0)²), darüber nicht (Basis 100: {mp.nstr(produkte['100'][9], 4)}); "
                        f"Re T ≥ 1 auf der rechten Halbebene (vier Basen); ∫_{{−1}}^0 T_b = 1/2 bei b = {mp.nstr(b_halb, 15)}; "
                        f"Carleman-Eigenwerte → λ^k (1.3: Abstand {mp.nstr(carl['1.3'][32], 2)} bei N = 32).")
lauf.schreibe()
