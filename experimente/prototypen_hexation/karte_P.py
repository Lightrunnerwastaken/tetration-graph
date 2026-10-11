"""Phasenbild der Hauptzweig-Pentation P_2 auf [−3, 3] × [−0.6, 4.2]: Schnitte als Farbsprünge, dazu Fixpunkte p, q und der Startweg H₀, P(H₀)."""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from blatt import BlattPentation
mp.mp.dps = 30
t0 = time.time()
BP = BlattPentation("2")
xs = np.arange(-3, 3.0001, 0.05)
ys = np.arange(-0.6, 4.2001, 0.05)
pkt = [mp.mpc(x, y) for y in ys for x in xs]
werte = BP.Pe.P(pkt)
Z = np.array([complex(w) for w in werte]).reshape(len(ys), len(xs))
print(f"P auf {len(pkt)} Punkten [{time.time() - t0:.0f}s]", flush=True)
np.save(HIER + "/prototypen_hexation/karte_P.npy", Z)
q = complex(-1.9929484359323223, 1.4451284966392055)
p = complex(1.7937303432391166, 1.5736302435118266)
st = json.load(open(HIER + "/prototypen_hexation/start_verfolgt_12.json"))
H = np.array([complex(a, b) for a, b in st["H"]])
fig, ax = plt.subplots(1, 2, figsize=(16, 7))
for a, f, titel in ((ax[0], np.angle(Z - q), "arg(P(w) − q)"), (ax[1], np.log10(np.abs(Z) + 1e-300).clip(-2, 4), "log10 |P(w)|")):
    im = a.imshow(f, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="hsv" if a is ax[0] else "viridis", aspect="equal")
    a.plot(H.real, H.imag, "k-", lw=1.5, label="Startweg H₀")
    a.plot([q.real], [q.imag], "wo", ms=7, mec="k"); a.plot([p.real], [p.imag], "ws", ms=7, mec="k")
    a.axhline(float(mp.pi / mp.re(BP.Pe.lnlam)), color="k", ls=":", lw=0.8)
    a.set_title(titel); plt.colorbar(im, ax=a, shrink=0.7)
plt.tight_layout()
plt.savefig(HIER + "/prototypen_hexation/karte_P.png", dpi=90)
print("Bild gespeichert", flush=True)
