"""exp.randgesetz — Sattel-Knoten-Durchgang nahe η und das Randgesetz von m(b).

(1) Reine Iteration (unabhängig von Kneser): N(b) = Schritte von t = 1 bis t > e + 1 unter t ↦ b^t.
    Vorhersage br.satz.sattel_knoten: N(b) = A·(b−η)^{−1/2} + O(1), A = π√(2η/e) = 3.23907…
(2) Engine-Daten m(b) = µ_{e,b} aus gp-tetration (mu_hub.json): Anpassung m = −A(b−η)^{−1/2} + B + C ln(b−η).
"""
import mpmath as mp

from gemeinsam import Lauf, gp_json, plt

mp.mp.dps = 40
lauf = Lauf("exp.randgesetz", __file__)
eta = mp.e ** (1 / mp.e)
A_th = mp.pi * mp.sqrt(2 * eta / mp.e)


def schritte(b):
    t, n = mp.mpf(1), 0
    while t <= mp.e + 1:
        t = mp.power(b, t)
        n += 1
    return n


abst = [mp.mpf(10) ** (-k / mp.mpf(2)) for k in range(4, 15)]   # b − η von 1e-2 bis 1e-7
N = [schritte(eta + d) for d in abst]
x = [d ** mp.mpf(-0.5) for d in abst]
# lineare Anpassung N = A x + B
n_ = len(x)
mx, my = sum(x) / n_, sum(N) / n_
A_fit = sum((xi - mx) * (yi - my) for xi, yi in zip(x, N)) / sum((xi - mx) ** 2 for xi in x)
B_fit = my - A_fit * mx
lauf.pruefe("Iterationszahl: Steigung A trifft π√(2η/e) auf 0.1 %", abs(A_fit / A_th - 1) < 1e-3,
            f"A = {mp.nstr(A_fit, 7)} vs {mp.nstr(A_th, 7)}")
rest = [Ni - A_th * xi for Ni, xi in zip(N, x)]
lauf.pruefe("Rest N − A·(b−η)^{−1/2} bleibt beschränkt (Spannweite < 3)", max(rest) - min(rest) < 3,
            [mp.nstr(r, 3) for r in rest])

hub = gp_json("mu_hub.json")["values"]
bs = sorted((mp.mpf(k), mp.mpf(v)) for k, v in hub.items() if mp.mpf(k) <= 1.8)
X = mp.matrix([[-(b - eta) ** mp.mpf(-0.5), 1, mp.log(b - eta)] for b, _ in bs])
Y = mp.matrix([m for _, m in bs])
koef = mp.lu_solve(X.T * X, X.T * Y)
A_m = koef[0]
lauf.pruefe("Engine-Hub m(b): angepasster Leitkoeffizient trifft π√(2η/e) auf 5 %", abs(A_m / A_th - 1) < 0.05,
            f"A = {mp.nstr(A_m, 6)} (B = {mp.nstr(koef[1], 4)}, C = {mp.nstr(koef[2], 4)})")
X2 = mp.matrix([[-(b - eta) ** mp.mpf(-0.5), 1] for b, _ in bs])
k2 = mp.lu_solve(X2.T * X2, X2.T * Y)
lauf.daten = {"N": dict((mp.nstr(d, 3), n) for d, n in zip(abst, N)), "A_iteration": A_fit, "B_iteration": B_fit,
              "A_hub_mit_log": A_m, "A_hub_ohne_log": k2[0], "A_theorie": A_th}

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot([float(v) for v in x], N, "o", label="Iterationszahl N(b)")
ax[0].plot([float(v) for v in x], [float(A_th * v + B_fit) for v in x], "--", label="π√(2η/e)·(b−η)^{−1/2} + B")
ax[0].set_xlabel("(b − η)^{−1/2}")
ax[0].set_ylabel("Schritte bis t > e + 1")
ax[0].set_title("Sattel-Knoten-Durchgang (reine Iteration)")
ax[0].legend(fontsize=8)
ax[1].plot([float((b - eta) ** -0.5) for b, _ in bs], [float(m) for _, m in bs], "o", label="m(b) aus gp-tetration")
gx = [mp.mpf(1.47) + i * mp.mpf("0.005") for i in range(70)]
ax[1].plot([float((g - eta) ** -0.5) for g in gx],
           [float(-koef[0] * (g - eta) ** mp.mpf(-0.5) + koef[1] + koef[2] * mp.log(g - eta)) for g in gx], "--",
           label="Anpassung −A x + B + C ln(b−η)")
ax[1].set_xlabel("(b − η)^{−1/2}")
ax[1].set_ylabel("m(b) = µ_{e,b}")
ax[1].set_title("Randgesetz der Hub-Spalte")
ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Iteration: A = {mp.nstr(A_fit, 6)} vs π√(2η/e) = {mp.nstr(A_th, 6)}. "
                        f"Hub-Daten: A = {mp.nstr(A_m, 5)} mit ln-Term, {mp.nstr(k2[0], 5)} ohne.")
lauf.schreibe()
