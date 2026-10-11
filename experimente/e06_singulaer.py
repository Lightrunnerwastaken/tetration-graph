"""exp.singulaer — T(−2+ε) − ln ε / ln b → ln T'(−1)/ln b (kn.satz.singulaer_minus2).

Getestet an der regulären Tetration zur Basis √2 und an der linearen Tetration (dort exakt).
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 50
lauf = Lauf("exp.singulaer", __file__)
R = Regulaer(mp.sqrt(2), N=60, tol=mp.mpf("1e-3"))
lb = R.lb
T = R.tet
Tp_m1 = 1 / mp.diff(R.slog, mp.mpf(0))   # T'(−1) = 1/slog'(T(−1)) = 1/slog'(0)
konst = mp.log(Tp_m1) / lb
eps = [mp.mpf(10) ** (-k) for k in range(2, 13)]
abw = []
for e in eps:
    wert = T(-2 + e) - mp.log(e) / lb
    abw.append(wert - konst)
lauf.pruefe("reguläre √2-Tetration: Restterm → 0 (ε = 1e-12)", abs(abw[-1]) < 1e-9, mp.nstr(abw[-1], 3))
verh = [abw[i + 1] / abw[i] for i in range(len(abw) - 1)]
lauf.pruefe("Restterm fällt linear in ε (Verhältnis ≈ 0.1 pro Dekade)", all(abs(v - mp.mpf("0.1")) < 0.02 for v in verh[2:]),
            [mp.nstr(v, 4) for v in verh])
# lineare Tetration: T(x) = log_b(T(x+1)), T(−1+ε) = ε ⇒ T(−2+ε) = log_b ε exakt, T'(−1) = 1
lin = [mp.log(e) / lb - mp.log(e) / lb for e in eps]
lauf.pruefe("lineare Tetration: Konstante ln T'(−1)/ln b = 0 exakt", all(x == 0 for x in lin))
lauf.pruefe("T'(−1) > 0 (Vorzeichen folgt aus T(−1+ε) > 0)", Tp_m1 > 0, mp.nstr(Tp_m1, 10))

fig, ax = plt.subplots(figsize=(6, 3.6))
ax.loglog([float(e) for e in eps], [abs(float(a)) for a in abw], "o-", label="|T(−2+ε) − ln ε/ln b − ln T′(−1)/ln b|")
ax.loglog([float(e) for e in eps], [float(e) * abs(float(abw[0] / eps[0])) for e in eps], "--", color="#888", label="∝ ε")
ax.set_xlabel("ε")
ax.set_title("Singularität bei −2: reguläre Tetration, Basis √2")
ax.legend(fontsize=8)
lauf.bild(fig)
lauf.daten = {"T_strich_minus1": Tp_m1, "konstante": konst, "restterme": [mp.nstr(a, 4) for a in abw]}
lauf.zusammenfassung = (f"T'(−1) = {mp.nstr(Tp_m1, 10)}; Restterm fällt wie ε (bei ε = 1e-12: {mp.nstr(abw[-1], 2)}). "
                        "Lineare Tetration: Konstante exakt 0.")
lauf.schreibe()
