"""exp.polynom_team — Nachrechnung der Kernzahlen des Forschungsteams (forschung/polynom_halbexp, 10.10.2026) aus dessen Daten.

 Daten (Strang A, gemeinsam/): Π_N f_K = J_N(f_K)⁻¹T_N[eˣ + (f_K'∘f_K)·f_K] (lineare J-Projektion) für N = 100 … 1000 (60 Stellen),
 echte Lösungen p_N für N = 100, 200, …, 1000, Knesers Taylor-Koeffizienten bis k = 1000 und die Singularitätsdaten an L.
 Geprüft wird nur, was sich aus diesen Zahlen direkt ablesen lässt:
   1. max_{|z|=|L|} |Π_N f_K| bleibt 2.232651 (Voraussetzung von Satz B-lin mit a = 0);
   2. Randfehler |Π_N f_K(L) − L| fällt wie N^{−4.8…−4.9};
   3. Innenfehler |Π_N f_K(1) − T_e(1/2)| fällt (Hüllkurve) wie N^{−7.9} — das Amplitudengesetz |c_N| ~ N^{1−2 Re γ};
   4. echte p_N liegen quadratisch nah an Π_N f_K: max|p_N − Π_N f_K| auf |z| = |L| ≤ 10·N·|Π_N f_K(L) − L|²;
   5. Knesers Koeffizienten folgen dem Darboux-Leitterm a_k ≈ 2 Re[D₁ L^{−k} Γ(k−1−γ)/(Γ(−1−γ) k!)] bis auf O(1/k).
"""
import json
from pathlib import Path

import mpmath as mp

from gemeinsam import Lauf, kneser, plt

mp.mp.dps = 60
lauf = Lauf("exp.polynom_team", __file__)
G = Path(__file__).resolve().parent.parent / "forschung" / "polynom_halbexp" / "gemeinsam"
lin = json.loads((G / "lin_PiN_fK_bis_1000.json").read_text(encoding="utf-8"))
echt = json.loads((G / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"]
sing = json.loads((G / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
koeff = json.loads((G / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))

L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
T_halb = mp.re(kneser("sexp", "exp(1)", [mp.mpf("0.5")])[0])


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


Ns = sorted(int(n) for n in lin)
Pi = {N: [mp.mpf(c) for c in lin[str(N)]["b"]] for N in Ns}
kreis = [abs(L) * mp.expj(2 * mp.pi * j / 720) for j in range(720)]
stuetz = [N for N in Ns if N % 100 == 0]
maxima = {N: max(abs(horner(Pi[N], z)) for z in kreis) for N in stuetz}
lauf.pruefe("max_{|z|=|L|} |Π_N f_K| = 2.232651 für N = 100 … 1000 (Satz B-lin braucht nur Wachstum < |z−L|^−3.447)",
            all(abs(v - mp.mpf("2.232651")) < mp.mpf("2e-6") for v in maxima.values()), {N: mp.nstr(v, 8) for N, v in maxima.items()})

rand = {N: abs(horner(Pi[N], L) - L) for N in Ns}
x = [mp.log(N) for N in Ns if N >= 200]
y = [mp.log(rand[N]) for N in Ns if N >= 200]
n = len(x)
steig = (n * mp.fsum(a * b for a, b in zip(x, y)) - mp.fsum(x) * mp.fsum(y)) / (n * mp.fsum(a * a for a in x) - mp.fsum(x) ** 2)
lauf.pruefe("Randfehler |Π_N f_K(L) − L| ~ N^p mit p ≈ −4.8 (Ecke von W mit Aussenwinkel 0.883π)", mp.mpf("-5.1") < steig < mp.mpf("-4.6"),
            {"p": mp.nstr(steig, 4), "N=200": mp.nstr(rand[200], 3), "N=1000": mp.nstr(rand[1000], 3)})

ak = [mp.mpf(s) for s in koeff["a"]]
kreis1 = [mp.expj(2 * mp.pi * j / 180) for j in range(180)]
fK1 = [horner(ak, z) for z in kreis1]
innen = {N: max(abs(horner(Pi[N], z) - w) for z, w in zip(kreis1, fK1)) for N in stuetz}
xs, ys = [mp.log(N) for N in stuetz if N >= 200], [mp.log(innen[N]) for N in stuetz if N >= 200]
m = len(xs)
steig_i = (m * mp.fsum(a * b for a, b in zip(xs, ys)) - mp.fsum(xs) * mp.fsum(ys)) / (m * mp.fsum(a * a for a in xs) - mp.fsum(xs) ** 2)
soll = 1 - 2 * mp.re(2j * mp.pi / L)
lauf.pruefe("Innenfehler max_{|z|=1} |Π_N f_K − f_K| ~ N^p mit p ≈ 1 − 2 Re(2πi/L) = −7.894 (Amplitudengesetz des Störmodus)",
            abs(steig_i - soll) < mp.mpf("0.4"), {"p": mp.nstr(steig_i, 4), "soll": mp.nstr(soll, 5), "N=1000": mp.nstr(innen[1000], 3)})

quad = {}
for N in sorted(int(k) for k in echt):
    p = [mp.mpf(c) for c in echt[str(N)]]
    quad[N] = (max(abs(horner(p, z) - horner(Pi[N], z)) for z in kreis[::6]), N * rand[N] ** 2)
lauf.pruefe("echte p_N (Newton, quadratisch) liegen bei Π_N f_K: max|p_N − Π_N f_K| ≤ 10·N·|E_N(L)|² auf |z| = |L| (N = 100 … 1000)",
            all(a <= 10 * b for a, b in quad.values()), {N: [mp.nstr(a, 3), mp.nstr(b, 3)] for N, (a, b) in quad.items() if N % 300 == 100})

D1 = mp.mpc(*[mp.mpf(s) for s in sing["D1_form_1_minus_z_durch_L"]])
ga = 2j * mp.pi / L
rel = {}
for k in (100, 200, 300, 400, 600):
    leit = 2 * D1 * L ** (-k) * mp.gamma(k - 1 - ga) / (mp.gamma(-1 - ga) * mp.factorial(k))
    rel[k] = abs(mp.re(leit / 2 * 2) - ak[k]) / abs(leit)          # Abstand relativ zur Leitterm-Amplitude (a_k schwingt durch 0)
lauf.pruefe("Knesers Koeffizienten = Darboux-Leitterm der Singularität (1 − z/L)^{1+γ}: Abstand ≤ 12/k der Leitterm-Amplitude (nächste Korrektur (1−z/L)^{2+γ} ist O(|1+γ|/k))",
            all(rel[k] < mp.mpf(12) / k for k in rel), {k: mp.nstr(v, 3) for k, v in rel.items()})

lauf.daten = {"max_kreis": {N: mp.nstr(v, 10) for N, v in maxima.items()}, "rand_exponent": mp.nstr(steig, 5),
              "innen_exponent": mp.nstr(steig_i, 5), "darboux_rel": {k: mp.nstr(v, 4) for k, v in rel.items()}}
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].loglog([float(N) for N in Ns], [float(rand[N]) for N in Ns], label="|Π_N f_K(L) − L|")
ax[0].loglog([float(N) for N in stuetz], [float(innen[N]) for N in stuetz], "o-", ms=3, label="max_{|z|=1}|Π_N f_K − f_K|")
ax[0].set_xlabel("N"); ax[0].legend(fontsize=8); ax[0].set_title("Rand- und Innenfehler der linearen Projektion")
ks = list(range(20, 1001, 5))
ax[1].semilogy(ks, [float(abs(ak[k]) * abs(L) ** k) for k in ks], lw=0.8, label="|a_k|·|L|^k (Kneser)")
ax[1].semilogy(ks, [71.57 * k ** -6.447 for k in ks], "k--", lw=0.8, label="71.57·k^−6.447")
ax[1].set_xlabel("k"); ax[1].legend(fontsize=8); ax[1].set_title("Knesers Koeffizienten und Singularität an L")
lauf.bild(fig)
lauf.zusammenfassung = (f"Team-Daten bestätigt: max|Π_N f_K| = 2.232651 bis N = 1000, Randfehler ~N^{mp.nstr(steig, 3)}, Innenfehler ~N^{mp.nstr(steig_i, 3)} "
                        "(Soll −7.894), echte p_N quadratisch nah, Kneser-Koeffizienten = Darboux-Leitterm an L.")
lauf.schreibe()
