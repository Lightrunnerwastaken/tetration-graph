"""Reelle Menge von F_q rechts der alten Karte (Re ζ bis 7) zwischen C und C + Π — Suche nach Log-Punkten, deren reeller Strahl schliesst.
 F_q(ζ) = pent^n(q + Ψ(μ^(ζ−n))) mit Hauptzweig-pent aus dem Fork (nur gelesen; Ψ = Umkehrreihe von ksig, zur Laufzeit gebildet).
 Erst Abgleich mit karte_reell.npz (dort mit Blattverfolgung von oben), dann das Gitter. Aufruf: karte_rechts.py [dx]
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fatou_backend import FatouGP, find_default_gp_exe

HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente/prototypen_hexation/"
dx = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
t0 = time.time()
ZETA_C = complex(0.4031293, -1.7097883)
PI_ = 2j * np.pi / np.log(complex(1.368932747339784991147632, 2.488144667623474107366636))
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)


def fq(zs, grenze=4):
    ex = [f"iferr(if(1, kv=alarm({grenze}, kFq({z.real:.15f}+I*{z.imag:.15f})); if(type(kv)==\"t_ERROR\", 10^9, kv)),E,10^9)" for z in zs]
    v = np.array([complex(x) for x in gp.eval_batch("2", ex)])
    v[np.abs(v) > 1e8] = np.nan
    return v


try:
    gp.eval_batch("2", ["hexinit(sexp(1),0,1)",
                        "if(1, kpsi = Pol(serreverse(ksig + O(x^(poldegree(ksig)+1)))); kfrs = abs(subst(ksig, x, krs))/2; 0)",
                        "if(1, kFq = ((z) -> if(1, my(s = exp(z*klog), n = 0, y); while(abs(s) > kfrs, s = s/klam; n++); y = kp + subst(kpsi, x, s); for(i = 1, n, y = pent(y)); y)); 0)"])
    # Abgleich mit der verfolgten Karte
    alt = np.load(HIER + "karte_reell.npz")
    pr = [(x, y) for x in (-2.0, -1.0, 0.0, 1.0) for y in (1.0, 0.0, -1.0, -1.4, -2.0, -2.4)]
    ref = []
    for x, y in pr:
        i = int(np.argmin(np.abs(alt["ys"] - y))); j = int(np.argmin(np.abs(alt["xs"] - x)))
        ref.append(alt["Z"][i, j])
    v = fq([complex(x, y) for x, y in pr])
    print(f"Abgleich mit karte_reell (verfolgt) [{time.time() - t0:.0f}s]:", flush=True)
    for (x, y), a, b in zip(pr, ref, v):
        print(f"   ζ = {x:+.1f}{y:+.1f}i: verfolgt {a:.8f}, Hauptkette {b:.8f}, Abstand {abs(a - b):.1e}", flush=True)
    xs = np.arange(-1.0, 7.0 + 1e-9, dx)
    ys = np.arange(-1.8, 1.3 + 1e-9, dx / 2)
    Z = np.full((len(ys), len(xs)), np.nan + 0j)
    for i, y in enumerate(ys):
        Z[i] = fq([complex(x, y) for x in xs])
        if i % 10 == 0:
            print(f"[{time.time() - t0:5.0f}s] Zeile {i + 1}/{len(ys)} (y = {y:+.2f}), {np.isnan(Z[i]).sum()} Ausfälle", flush=True)
finally:
    gp.close()
np.savez(HIER + "karte_rechts.npz", xs=xs, ys=ys, Z=Z)
fig, ax = plt.subplots(1, 1, figsize=(15, 6.5))
ax.imshow(np.angle(Z), origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="hsv", aspect="equal")
ax.contour(xs, ys, Z.imag, levels=[0], colors="k", linewidths=1.0)
for k in (0, 1, 2):
    for m in range(-2, 9):
        z = ZETA_C - 2 + m + k * PI_
        if xs[0] <= z.real <= xs[-1] and ys[0] <= z.imag <= ys[-1]:
            ax.plot(z.real, z.imag, "wo" if m == 0 else "w^", ms=6, mec="k")
ax.set_title("arg F_q (Hauptkette), schwarz: Im F_q = 0; Kreis: Log-Punkt, Dreieck: Ecke (Zeilen k = 0, 1)")
plt.tight_layout()
plt.savefig(HIER + "karte_rechts.png", dpi=85)
print(f"fertig [{time.time() - t0:.0f}s]", flush=True)
