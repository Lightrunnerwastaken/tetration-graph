"""exp.young — Flächenidentität ∫(H−id) = −∫(H⁻¹−id) für Lifts von Kreishomöomorphismen.

(1) 200 zufällige glatte Lifts: Identität bis auf Quadraturfehler (Gegenbeispielsuche).
(2) Nicht-glatter Lift (stückweise linear): Identität gilt trotzdem (Glattheit unnötig).
(3) Anwendung auf echte Phasen Φ_{e,b} aus gp-tetration: ∫(H⁻¹−id) = −µ_{e,b} = µ_{b,e}.
"""
import random

import numpy as np

from gemeinsam import Lauf, gp_json, plt

lauf = Lauf("exp.young", __file__)
G = 4096
th = (np.arange(G) + 0.5) / G


def inverse_mittel(H_werte):
    """∫_0^1 (H⁻¹(y) − y) dy über die Umkehrung auf einem feinen Gitter."""
    # H auf [0,1) mit H(x+1) = H(x)+1; Umkehrung durch Interpolation auf erweitertem Gitter
    verschiebung = range(-6, 7)          # deckt auch grosse Mittel |µ| ≤ 5 ab
    xs = np.concatenate([th + s for s in verschiebung])
    hs = np.concatenate([H_werte + s for s in verschiebung])
    ys = th
    return float(np.mean(np.interp(ys, hs, xs) - ys))


random.seed(3)
worst = 0.0
for _ in range(200):
    mu = random.uniform(-0.5, 0.5)
    koef = [(random.uniform(-1, 1), random.uniform(0, 1)) for _ in range(4)]
    s = sum(abs(a) * 2 * np.pi * (k + 1) for k, (a, _) in enumerate(koef))
    skal = random.uniform(0.1, 0.9) / s          # sorgt für H' > 0
    H = th + mu + skal * sum(a * np.sin(2 * np.pi * (k + 1) * (th + ph)) for k, (a, ph) in enumerate(koef))
    worst = max(worst, abs(np.mean(H - th) + inverse_mittel(H)))
lauf.pruefe("200 zufällige glatte Lifts: |∫(H−id) + ∫(H⁻¹−id)| < 1e-6", worst < 1e-6, f"{worst:.2e}")

# nicht glatt: stückweise linearer Lift
H = np.where(th < 0.3, th * 0.5, 0.15 + (th - 0.3) * (0.85 / 0.7)) + 0.2
lauf.pruefe("stückweise linearer Lift (nicht C¹): Identität gilt", abs(np.mean(H - th) + inverse_mittel(H)) < 1e-6,
            f"{abs(np.mean(H - th) + inverse_mittel(H)):.2e}")

modes = gp_json("phi_modes.json")["bases"]
tab = []
for b, dat in modes.items():
    mu = float(dat["mu"])
    Phi = mu + sum(float(m["re"]) * np.cos(2 * np.pi * (k + 1) * th) - float(m["im"]) * np.sin(2 * np.pi * (k + 1) * th)
                   for k, m in enumerate(dat["modes"]))
    H = th + Phi
    inv = inverse_mittel(H)
    tab.append({"b": b, "mu_e_b": mu, "inverses_mittel": inv, "summe": mu + inv})
abw = max(abs(t["summe"]) for t in tab)
lauf.pruefe(f"Phasen Φ_(e,b) aus gp-tetration ({len(tab)} Basen): ∫(H⁻¹−id) = −µ_(e,b)", abw < 1e-6, f"{abw:.2e}")
lauf.daten = {"zufall_max_abw": worst, "phasen": tab}

fig, ax = plt.subplots(figsize=(5.5, 3.6))
ax.bar([t["b"] for t in tab], [abs(t["summe"]) + 1e-18 for t in tab], color="#3d6fb6")
ax.set_yscale("log")
ax.set_xlabel("Basis b")
ax.set_ylabel("|µ_(e,b) + ∫(H⁻¹ − id)|")
ax.set_title("Flächenidentität an den echten Phasen (Quadraturfehler)")
lauf.bild(fig)
lauf.zusammenfassung = (f"Identität für 200 Zufallslifts (max. Abw. {worst:.1e}), einen nicht-glatten Lift und {len(tab)} echte "
                        f"Phasen (max. {abw:.1e}) bestätigt: Paper II Thm. B braucht nur Monotonie + Inverse.")
lauf.schreibe()
