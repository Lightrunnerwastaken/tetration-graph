"""exp.ableitungen_vorzeichen — Vorzeichen der Ableitungen der Kneser-Tetration T_e (Forum tid 882, mike3 / Levenstein 2014).

 Vermutung (mike3, 09.06.2014): alle ungeraden Ableitungen von T_e sind auf (−2, ∞) positiv (gleichwertig: alle geraden wachsen).
 Levenstein (14.06.2014): hohe Ableitungen folgen log(x+2) (nächste Singularität −2), daher sind die geraden Ableitungen hoher Ordnung
 überall negativ, z. B. bei z = 0.5 zuerst die 406. — die geraden wechseln also ihr Vorzeichen.
 Hier: Taylor-Koeffizienten a_k(x0) = T^{(k)}(x0)/k! für k ≤ 26 an 20 Stellen x0 ∈ [−1.8, 2]; Vorzeichenwechsel der geraden Ordnungen;
 Vergleich mit log(x+2): a_k·k·(x0+2)^k·(−1)^{k+1} → 1.
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from stufen import taylor

mp.mp.dps = 40
lauf = Lauf("exp.ableitungen_vorzeichen", __file__)
K = 26


def koeff(x0):
    r = min(mp.mpf("0.8") * (x0 + 2), mp.mpf(1))
    return taylor("exp(1)", x0, r, N=96, K=K)


gitter = [mp.mpf(-18 + 2 * i) / 10 for i in range(20)]           # −1.8 … 2.0
A = {x: koeff(x) for x in gitter}
ungerade_ok = all(A[x][k] > 0 for x in gitter for k in range(1, K + 1, 2))
lauf.pruefe("alle ungeraden Ableitungen bis Ordnung 25 positiv an 20 Stellen in [−1.8, 2]", ungerade_ok,
            {mp.nstr(x, 2): mp.nstr(min(A[x][k] for k in range(1, K + 1, 2)), 3) for x in gitter[::4]})

# gerade Ordnungen: Vorzeichenfolge entlang des Gitters
wechsel = {}
for k in range(2, K + 1, 2):
    vz = [A[x][k] > 0 for x in gitter]
    n = sum(1 for i in range(len(vz) - 1) if vz[i] != vz[i + 1])
    wechsel[k] = (n, vz[0], vz[-1])
lauf.pruefe("jede gerade Ableitung (2 … 26) wechselt auf dem Gitter genau einmal, von − nach +",
            all(n == 1 and not links and rechts for n, links, rechts in wechsel.values()), {k: v[0] for k, v in wechsel.items()})


# Wechselstellen x_{2n} genauer (Sekanten in x0)
def nullstelle(k):
    i = next(i for i in range(len(gitter) - 1) if (A[gitter[i]][k] > 0) != (A[gitter[i + 1]][k] > 0))
    x0, x1 = gitter[i], gitter[i + 1]
    f0, f1 = A[x0][k], A[x1][k]
    for _ in range(40):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
        f1 = koeff(x1)[k]
        if abs(x1 - x0) < mp.mpf(10) ** -20:
            break
    return x1


xs_w = {k: nullstelle(k) for k in (2, 4, 6, 10, 16, 22, 26)}
lauf.pruefe("x₂ (Wendepunkt) = −0.5218… wie gr.beob.wendepunkte", abs(xs_w[2] - mp.mpf("-0.5218")) < mp.mpf("5e-4"), mp.nstr(xs_w[2], 15))
lauf.pruefe("die Wechselstellen x_{2n} wandern nach rechts", all(xs_w[a] < xs_w[b] for a, b in zip(list(xs_w)[:-1], list(xs_w)[1:])),
            {k: mp.nstr(v, 10) for k, v in xs_w.items()})
# log(x+2)-Gesetz bei x0 = −1: a_k ≈ (−1)^{k+1} / (k (x0+2)^k)
a_m1 = A[mp.mpf(-1)] if mp.mpf(-1) in A else koeff(mp.mpf(-1))
quot = {k: a_m1[k] * k * (-1) ** (k + 1) for k in (6, 12, 18, 24)}
lauf.pruefe("hohe Ableitungen folgen log(x+2): a_k·k·(−1)^{k+1} → 1 bei x0 = −1 (Abstand wie 2^{−k})",
            abs(quot[24] - 1) < mp.mpf(10) ** -5 and abs(quot[24] - 1) < abs(quot[12] - 1) < abs(quot[6] - 1), {k: mp.nstr(v, 10) for k, v in quot.items()})
# Levensteins Faustregel: Ordnung des ersten negativen Koeffizienten bei 0.5 ist 406 — hier nur: bis 26 noch alle gerade positiv
lauf.pruefe("bei x0 = 0.5 sind alle Ableitungen bis Ordnung 26 positiv (Levenstein: erste negative ist die 406.)",
            all(koeff(mp.mpf("0.5"))[k] > 0 for k in range(1, K + 1)))

lauf.daten = {"wechselstellen": {k: mp.nstr(v, 20) for k, v in xs_w.items()},
              "log_quotient_bei_minus1": {k: mp.nstr(v, 15) for k, v in quot.items()},
              "min_ungerade": {mp.nstr(x, 3): mp.nstr(min(A[x][k] for k in range(1, K + 1, 2)), 5) for x in gitter}}

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for k, f in ((2, "#9fb7d9"), (6, "#6f93c6"), (10, "#3d6fb6"), (16, "#244f8f"), (22, "#b5523b"), (26, "#8f3a26")):
    ys = [A[x][k] * (x + 2) ** k for x in gitter]                 # skaliert, damit alle Ordnungen in ein Bild passen
    ax[0].plot([float(x) for x in gitter], [float(mp.sign(y) * mp.log10(1 + abs(y) * 1e3)) for y in ys], "o-", ms=2.5, color=f, label=f"Ordnung {k}")
ax[0].axhline(0, color="#8a93a3", lw=0.6)
ax[0].set_xlabel("x"); ax[0].set_ylabel("sign · log₁₀(1 + 10³·|a_k|(x+2)^k)"); ax[0].set_title("Gerade Ableitungen: ein Vorzeichenwechsel"); ax[0].legend(fontsize=7)
ax[1].plot(list(xs_w), [float(v) for v in xs_w.values()], "o-", color="#b5523b")
ax[1].set_xlabel("Ordnung 2n"); ax[1].set_ylabel("Wechselstelle x_{2n}"); ax[1].set_title("Die Wechselstellen wandern nach rechts")
lauf.bild(fig)
lauf.zusammenfassung = (f"Ungerade Ableitungen bis 25 überall positiv (20 Stellen); gerade wechseln genau einmal: x₂ = {mp.nstr(xs_w[2], 8)}, "
                        f"x₁₀ = {mp.nstr(xs_w[10], 6)}, x₂₂ = {mp.nstr(xs_w[22], 6)}, x₂₆ = {mp.nstr(xs_w[26], 6)}; log(x+2)-Gesetz bestätigt.")
lauf.schreibe()
