"""exp.nichteindeutigkeit — zwei reell-analytische Tetrationen zur Basis √2 mit gleichem T(0)=1.

F = reguläre Tetration (Koenigs am Fixpunkt 2). T̃ = F∘(id+θ), θ(x) = c(1 − cos 2πx)/2,
c so gewählt, dass T̃(1/2) = r für vorgegebenes r. Geprüft: Gleichung, Normierung, Monotonie.
"""
import random

import mpmath as mp

from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 50
lauf = Lauf("exp.nichteindeutigkeit", __file__)
R = Regulaer(mp.sqrt(2), N=60, tol=mp.mpf("1e-3"))
b = R.b


def F(x):
    return R.tet(x)


def theta(x, c):
    return c * (1 - mp.cos(2 * mp.pi * x)) / 2


fall = []
for r in (mp.mpf("1.15"), mp.mpf("1.30")):
    c = R.slog(r) - mp.mpf("0.5")
    lauf.pruefe(f"r={r}: |c| < 1/π (id+θ streng wachsend)", abs(c) < 1 / mp.pi, mp.nstr(c, 8))

    def Tt(x, c=c):
        return F(x + theta(x, c))

    random.seed(1)
    res = max(abs(Tt(x + 1) - mp.power(b, Tt(x))) for x in [mp.mpf(random.uniform(-1.5, 2)) for _ in range(12)])
    lauf.pruefe(f"r={r}: T̃(x+1) = b^T̃(x) an 12 Zufallspunkten", res < mp.mpf(10) ** -40, mp.nstr(res, 3))
    lauf.pruefe(f"r={r}: T̃(0) = 1 und T̃(1/2) = r", abs(Tt(0) - 1) < 1e-40 and abs(Tt(mp.mpf("0.5")) - r) < 1e-40,
                mp.nstr(Tt(mp.mpf("0.5")), 12))
    ableit = min(mp.diff(Tt, x, h=mp.mpf("1e-15")) for x in [mp.mpf(i) / 20 for i in range(-30, 40)])
    lauf.pruefe(f"r={r}: T̃' > 0 auf [−1.5, 2]", ableit > 0, mp.nstr(ableit, 5))
    fall.append((r, c, Tt))

F_halb = F(mp.mpf("0.5"))
lauf.pruefe("die drei Lösungen unterscheiden sich bei 1/2", abs(F_halb - fall[0][0]) > 0.01 and abs(F_halb - fall[1][0]) > 0.01,
            mp.nstr(F_halb, 10))
xs = [mp.mpf(i) / 40 for i in range(-70, 121)]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot([float(x) for x in xs], [float(F(x)) for x in xs], label="reguläre Tetration F")
for r, c, Tt in fall:
    ax[0].plot([float(x) for x in xs], [float(Tt(x)) for x in xs], label=f"F∘(id+θ), T̃(1/2) = {r}")
    ax[1].plot([float(x) for x in xs], [float(Tt(x) - F(x)) for x in xs], label=f"T̃ − F (r = {r})")
ax[0].axhline(2, ls=":", color="#888")
ax[0].set_xlabel("x")
ax[0].set_title("Drei reell-analytische Tetrationen zur Basis √2")
ax[0].legend(fontsize=8)
ax[1].set_xlabel("x")
ax[1].set_title("Differenz zur regulären Lösung")
ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.daten = {"F(1/2)": F_halb, "faelle": [{"r": r, "c": c} for r, c, _ in fall]}
lauf.zusammenfassung = (f"Reguläre Tetration F(1/2) = {mp.nstr(F_halb, 12)}; zwei weitere reell-analytische, streng wachsende "
                        "Lösungen mit T̃(1/2) = 1.15 bzw. 1.30 erfüllen die Gleichung auf 1e-40.")
lauf.schreibe()
