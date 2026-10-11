"""exp.polynom_halbexp — Polynom-Wurzel von exp: wählt das Abschneiden die Kneser-Halbexponentialfunktion aus?

 Verfahren (Janis' Notiz, Okt. 2026): p_N(x) = Σ_{k≤N} a_k x^k mit p_N(p_N(x)) = e^x + O(x^{N+1}) — N+1 Gleichungen für N+1 Koeffizienten,
 Zweig stetig von p_1(x) = x + 1/2 aus fortgesetzt (Newton, Startwert: Lösung zu N−1 mit a_N = 0). Gerechnet mit arb (python-flint).
 Vergleich: Kneser f_K(z) = T_e(slog_e(z) + 1/2), Taylor-Koeffizienten um 0 aus Cauchy-Werten auf |z| = 0.9 (komplexes slog und sexp).
 Prüfgrössen: p_N(1) gegen T_e(1/2), Koeffizienten einzeln, und p_N(L) gegen L (L = Fixpunkt von exp, |L| = 1.3746 = Konvergenzradius von f_K).
"""
import mpmath as mp

from gemeinsam import Lauf, kneser, plt
from polynom_wurzel import horner, loesungen

mp.mp.dps = 45
lauf = Lauf("exp.polynom_halbexp", __file__)
NMAX = 100
loes = loesungen(NMAX)

T_halb, T_mhalb = [mp.re(v) for v in kneser("sexp", "exp(1)", [mp.mpf("0.5"), mp.mpf("-0.5")])]
wert1 = {N: mp.fsum(c) for N, c in loes.items()}
abst = {N: wert1[N] - T_halb for N in loes}
NOTIZ = {2: "1.637802440995079", 3: "1.646465116557290", 5: "1.646354339729278", 15: "1.646354233718347"}
lauf.pruefe("die Werte p_N(1) der Notiz (N = 2, 3, 5, 15) sind nachgerechnet",
            all(abs(wert1[N] - mp.mpf(v)) < mp.mpf("1e-15") for N, v in NOTIZ.items()), {N: mp.nstr(wert1[N], 16) for N in NOTIZ})
lauf.pruefe("p_N(1) → T_e(1/2): Abstand bei N = 100 unter 2·10⁻¹⁶", abs(abst[100]) < mp.mpf("2e-16"),
            {N: mp.nstr(abst[N], 3) for N in (5, 10, 20, 40, 60, 80, 100)})
rate = (abs(abst[100] + abst[99]) / abs(abst[70] + abst[69])) ** (mp.mpf(1) / 30)
lauf.pruefe("ab N ≈ 50 geometrisch: Abstand schrumpft etwa um den Faktor 0.93 pro Grad (gemittelt über gerade/ungerade N)",
            mp.mpf("0.9") < rate < mp.mpf("0.95"), mp.nstr(rate, 4))

# Kneser-Koeffizienten um 0
r, M, K = mp.mpf("0.9"), 256, 80
zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
fz = kneser("sexp", "exp(1)", [w + mp.mpf("0.5") for w in kneser("slog", "exp(1)", zs)])
kn = [mp.re(mp.fsum(fz[j] * mp.expjpi(-2 * k * mp.mpf(j) / M) for j in range(M)) / M) / r ** k for k in range(K + 1)]
lauf.pruefe("Kneser-Reihe um 0: a₀ = T_e(−1/2), Σ a_k = T_e(1/2) (81 Glieder)", abs(kn[0] - T_mhalb) < mp.mpf("1e-30")
            and abs(mp.fsum(kn) - T_halb) < mp.mpf("1e-18"), {"a0": mp.nstr(kn[0], 20), "Summe-T(1/2)": mp.nstr(mp.fsum(kn) - T_halb, 3)})
rel = {k: {N: abs(loes[N][k] - kn[k]) / abs(kn[k]) for N in (20, 50, 100)} for k in (0, 1, 5, 10, 20)}
lauf.pruefe("jeder Koeffizient a_k^(N) → a_k^Kneser (k = 0, 1, 5, 10, 20; N = 20, 50, 100)",
            all(v[100] < v[50] < v[20] for v in rel.values()), {k: [mp.nstr(x, 2) for x in v.values()] for k, v in rel.items()})

# Am Rand: p_N(L) → L = f_K(L)
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
rand ={N: abs(horner(loes[N], L) - L) for N in (20, 40, 60, 80, 100)}
lauf.pruefe("am Randpunkt L (Fixpunkt von exp, |L| = 1.3746): |p_N(L) − L| fällt (Kneser: f(L) = L)",
            all(rand[n2] < rand[n1] for n1, n2 in zip(list(rand)[:-1], list(rand)[1:])), {N: mp.nstr(v, 3) for N, v in rand.items()})

# Andere halbe e-Funktionen: T(φ⁻¹(φ(A)+½)) mit φ = id + ε·sin(2π·) — der Störterm ε·sin(2πA(z)) wächst zum Rand hin wie |z − L|^(−Re(2πi/L))
ts = [mp.mpf("0.9"), mp.mpf("0.99"), mp.mpf("0.999")]
A_rand = kneser("slog", "exp(1)", [t * L for t in ts])
stoer = [abs(mp.mpf("0.01") * mp.sin(2 * mp.pi * w)) for w in A_rand]
expo = [mp.log(stoer[i + 1] / stoer[i]) / mp.log(10) for i in range(2)]
soll = mp.re(2j * mp.pi / L)
lauf.pruefe("Störterm 0.01·sin(2πA(z)) anderer Lösungen wächst bei z = tL (t = 0.9, 0.99, 0.999) wie |z−L|^(−4.447), 4.447 = Re(2πi/L)",
            all(abs(x - soll) < mp.mpf("0.15") for x in expo), {"Störterm": [mp.nstr(s, 3) for s in stoer], "Exponent": [mp.nstr(x, 4) for x in expo],
                                                                 "Re(2πi/L)": mp.nstr(soll, 5)})

lauf.daten = {"p_N(1)": {N: mp.nstr(wert1[N], 30) for N in (1, 2, 3, 5, 10, 20, 50, 100)},
              "abstand": {N: mp.nstr(abst[N], 4) for N in range(1, NMAX + 1)},
              "koeff_N100": [mp.nstr(c, 25) for c in loes[100][:21]], "koeff_kneser": [mp.nstr(c, 25) for c in kn[:21]],
              "rate": mp.nstr(rate, 5), "rand_L": {N: mp.nstr(v, 4) for N, v in rand.items()}}
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ns = list(range(1, NMAX + 1))
ax[0].semilogy(ns, [max(float(abs(abst[N])), 1e-40) for N in ns], ".-", ms=3)
ax[0].set_xlabel("Grad N"); ax[0].set_title("|p_N(1) − T_e(1/2)|")
for k, v in rel.items():
    ax[1].semilogy([N for N in range(max(k, 1), NMAX + 1)], [float(abs(loes[N][k] - kn[k]) / abs(kn[k])) for N in range(max(k, 1), NMAX + 1)], label=f"k = {k}")
ax[1].set_xlabel("Grad N"); ax[1].set_title("relativer Abstand a_k^(N) zu Kneser"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Das Abschneiden p_N∘p_N = exp + O(x^(N+1)) trifft Kneser: |p_N(1) − T_e(1/2)| = {mp.nstr(abs(abst[100]), 2)} bei N = 100, "
                        f"Rate ≈ {mp.nstr(rate, 3)} pro Grad; jeder Koeffizient und der Randwert p_N(L) → L konvergieren.")
lauf.schreibe()
