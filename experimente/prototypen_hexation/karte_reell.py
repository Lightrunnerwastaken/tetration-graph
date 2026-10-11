"""Reelle Menge von F_q (Basis 2) mit Blattverfolgung: Spalten senkrecht von Im ζ = 3 (Hauptzweig dort richtig) nach unten.
 Frage zur Eindeutigkeit: Welche +1-invarianten Kurven mit reellem F_q gibt es ausser der Kneser-Kurve C (Ecken ζ_c − 2 + m)?
 Ausgabe: karte_reell.npz (Gitter, Werte) und karte_reell.png. Aufruf: karte_reell.py dx
"""
import sys, time
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from blatt import BlattPentation, verfolge
from hexq import HexQ
from stufen import psi_wert
mp.mp.dps = 30
t0 = time.time()
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
BP = BlattPentation("2")
PI_ = 2j * mp.pi / HQ.lnmu
zeta_c = mp.mpc("0.4031293", "-1.7097883")
dx = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
xs = np.arange(-2.6, 1.4 + 1e-9, dx)
ys = np.arange(-2.6, 1.6 + 1e-9, dx)
OBEN = 3.0


def F_weg(zetas, n):
    w = [psi_wert(HQ.q, HQ.psi, mp.exp((z - n) * HQ.lnmu)) for z in zetas]
    sp = [mp.mpf(0)] * len(zetas)
    for _ in range(n):
        w, s_ = BP.P_pfad(w)
        sp = [max(a, b) for a, b in zip(sp, s_)]
    return w, sp


Z = np.full((len(ys), len(xs)), np.nan + 0j)
for j, x in enumerate(xs):
    weg = [mp.mpc(x, OBEN - k * dx) for k in range(int(round((OBEN - ys[0]) / dx)) + 1)]
    n = max(0, max(int(mp.ceil(mp.re((mp.log(abs(mp.exp(z * HQ.lnmu))) - mp.log(HQ.rlok)) / mp.re(HQ.lnmu)))) for z in weg)) + 1
    try:
        p, w, sm = verfolge(weg, lambda zz: F_weg(zz, n), runden=7)
    except Exception as e:                                   # Spalte durch einen Verzweigungspunkt o. ä.
        print(f"  Spalte x = {x:.3f}: {type(e).__name__}", flush=True)
        continue
    werte = {round(float(mp.im(z)), 6): complex(v) for z, v in zip(p, w)}
    for i, y in enumerate(ys):
        v = werte.get(round(float(y), 6))
        if v is not None:
            Z[i, j] = v
    if j % 10 == 0:
        print(f"[{time.time() - t0:6.0f}s] Spalte {j + 1}/{len(xs)} (x = {x:.2f}), n = {n}, Sprung {mp.nstr(sm, 3)}", flush=True)
np.savez(HIER + "/prototypen_hexation/karte_reell.npz", xs=xs, ys=ys, Z=Z)
LT, PLT = -1.743909176132737, -0.771546447846277
fig, ax = plt.subplots(1, 2, figsize=(15, 7.5))
ax[0].imshow(np.angle(Z), origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="hsv", aspect="equal")
ax[0].contour(xs, ys, Z.imag, levels=[0], colors="k", linewidths=1.2)
ax[0].set_title("arg F_q, schwarz: Im F_q = 0 (reelle Menge)")
re = np.clip(Z.real, -4, 6)
re[np.abs(Z.imag) > 0.05 * (1 + np.abs(Z.real))] = np.nan
im1 = ax[1].imshow(re, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="viridis", aspect="equal")
ax[1].contour(xs, ys, Z.imag, levels=[0], colors="k", linewidths=0.8)
ax[1].set_title("Re F_q nahe der reellen Menge")
plt.colorbar(im1, ax=ax[1], shrink=0.7)
for a in ax:
    for k in (0, 1):
        for m in range(-1, 5):
            z = complex(zeta_c - 2 + m + k * PI_)
            if xs[0] <= z.real <= xs[-1] and ys[0] <= z.imag <= ys[-1]:
                a.plot(z.real, z.imag, "wo" if m == 0 else "w^", ms=6, mec="k")
plt.tight_layout()
plt.savefig(HIER + "/prototypen_hexation/karte_reell.png", dpi=90)
print(f"fertig [{time.time() - t0:.0f}s]", flush=True)
