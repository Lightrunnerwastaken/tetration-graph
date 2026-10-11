"""Reelle Menge von F_q UNTERHALB der Kneser-Kurve C (ein Streifen der Höhe Im Π): Zeilen waagrecht von Re ζ = −7 (dort F_q ≈ q,
 Hauptzweig richtig) nach rechts, mit Blattverfolgung. Schnitte von Singularitäten laufen so nach links; die Zeilen unter C kreuzen keine
 der Ecken von C. Frage: Gibt es dort Log-Punkte von F_q, deren reeller Strahl eine +1-invariante Kurve bildet?
 Aufruf: karte_unten.py dy
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
dy = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
DX = 0.05
LINKS = -7.0
xs = np.arange(-3.0, 1.6 + 1e-9, DX)
ys = np.arange(-4.6, -1.8 + 1e-9, dy)


def F_weg(zetas, n):
    w = [psi_wert(HQ.q, HQ.psi, mp.exp((z - n) * HQ.lnmu)) for z in zetas]
    sp = [mp.mpf(0)] * len(zetas)
    for _ in range(n):
        w, s_ = BP.P_pfad(w)
        sp = [max(a, b) for a, b in zip(sp, s_)]
    return w, sp


Z = np.full((len(ys), len(xs)), np.nan + 0j)
SPR = np.zeros(len(ys))
for i, y in enumerate(ys):
    weg = [mp.mpc(LINKS + k * DX, y) for k in range(int(round((xs[-1] - LINKS) / DX)) + 1)]
    n = max(0, max(int(mp.ceil(mp.re((mp.log(abs(mp.exp(z * HQ.lnmu))) - mp.log(HQ.rlok)) / mp.re(HQ.lnmu)))) for z in weg)) + 1
    try:
        p, w, sm = verfolge(weg, lambda zz: F_weg(zz, n), runden=7)
    except Exception as e:
        print(f"  Zeile y = {y:.2f}: {type(e).__name__} {str(e)[:60]}", flush=True)
        continue
    SPR[i] = float(sm)
    werte = {round(float(mp.re(z)), 6): complex(v) for z, v in zip(p, w)}
    for j, x in enumerate(xs):
        v = werte.get(round(float(x), 6))
        if v is not None:
            Z[i, j] = v
    print(f"[{time.time() - t0:6.0f}s] Zeile {i + 1}/{len(ys)} (y = {y:.2f}), n = {n}, {len(p)} Punkte, Sprung {mp.nstr(sm, 3)}", flush=True)
np.savez(HIER + "/prototypen_hexation/karte_unten.npz", xs=xs, ys=ys, Z=Z, spr=SPR)
fig, ax = plt.subplots(1, 2, figsize=(15, 7.5))
ax[0].imshow(np.angle(Z), origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="hsv", aspect="equal")
ax[0].contour(xs, ys, Z.imag, levels=[0], colors="k", linewidths=1.2)
ax[0].set_title("arg F_q unter C, schwarz: Im F_q = 0")
lg = np.log10(np.abs(Z) + 1e-300).clip(-2, 3)
im1 = ax[1].imshow(lg, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="viridis", aspect="equal")
ax[1].contour(xs, ys, Z.imag, levels=[0], colors="w", linewidths=0.8)
ax[1].set_title("log10 |F_q| (Log-Punkte: hell)")
plt.colorbar(im1, ax=ax[1], shrink=0.7)
for a in ax:
    for k in (0, -1):
        for m in range(-2, 5):
            z = complex(zeta_c - 2 + m + k * PI_)
            if xs[0] <= z.real <= xs[-1] and ys[0] - 0.3 <= z.imag <= ys[-1] + 0.3:
                a.plot(z.real, z.imag, "wo" if m == 0 else "w^", ms=6, mec="k")
plt.tight_layout()
plt.savefig(HIER + "/prototypen_hexation/karte_unten.png", dpi=90)
print(f"fertig [{time.time() - t0:.0f}s]", flush=True)
