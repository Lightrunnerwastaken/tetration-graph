"""Tabellen und Bild aus daten/*.json (keine neue Rechnung ausser den Fixpunkten L_b).
Aufruf: python zusammenfassung.py   ->  Tabellen auf stdout, Bild wellen_raenge.png"""
import json, math, glob
from pathlib import Path
import mpmath as mp
H = Path(__file__).resolve().parent
D = H / "daten"
mp.mp.dps = 40
def lade(pre, b, d, M):
    f = D / f"{pre}_{b}_{d}_M{M}.json"
    return json.loads(f.read_text(encoding="utf-8")) if f.is_file() else None
def wert(name):
    return {"exp1": mp.e}.get(name, None) or mp.mpf(eval(name.replace("_", "/")))
def ab(c):
    return math.hypot(float(c[0]), float(c[1]))
print("Rang 4 und 5 (E1):  b  d  M  ln(alpha)  mu4  mu5  |a1|_4  |a1|_5  Quotient  |a1|_5/|ln alpha|")
for f in sorted(D.glob("e1_*_M*.json")):
    r = json.loads(f.read_text(encoding="utf-8"))
    bb = {"exp(1)": mp.e}.get(r["b"]) or mp.mpf(eval(r["b"]))
    dd = {"exp(1)": mp.e}.get(r["d"]) or mp.mpf(eval(r["d"]))
    la = mp.log(mp.log(dd) / mp.log(bb))
    a4, a5 = ab(r["f4"][1]), ab(r["f5"][1])
    print(f"  {r['b']:>7} {r['d']:>6} {r['M']:>4} {float(la):+.5f}  {float(mp.mpf(r['mu4'])):+.12f} {float(mp.mpf(r['mu5'])):+.10f}  "
          f"{a4:.4e} {a5:.4e}  {a5/a4:6.2f}  {a5/abs(float(la)):.5f}")
print("Rang 6 (E2) und 7 (E3): b d M mu |a1| |a2| |a3| Spannweite")
for pre, k in (("e2", "6"), ("e3", "7")):
    for f in sorted(D.glob(f"{pre}_*_M*.json")):
        r = json.loads(f.read_text(encoding="utf-8"))
        fk = r["f" + k]
        print(f"  Rang {k}: {r['b']:>7} {r['d']:>6} {r['M']:>4} mu={float(mp.mpf(r['mu'+k])):+.10f} |a1|={ab(fk[1]):.4e} |a2|={ab(fk[2]):.4e} |a3|={ab(fk[3]):.4e} pp={float(mp.mpf(r['pp'+k])):.4e}")
# Bild: Wellen Phi_r - mu_r für (e,2), Rang 4..7
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for pre, k, M, farbe in (("e1", "4", 128, "#3d6fb6"), ("e1", "5", 128, "#b5523b"), ("e2", "6", 64, "#2f8f83"), ("e3", "7", 32, "#7a49c2")):
        r = lade(pre, "exp1", "2", M)
        v = [float(mp.mpf(x)) for x in r["phi" + k]]
        mu = sum(v) / len(v)
        xs = [j / len(v) for j in range(len(v) + 1)]
        ys = [(x - mu) * (50 if k == "4" else 1) for x in v + v[:1]]
        ax[0].plot(xs, ys, color=farbe, label=f"Rang {k}" + (" (×50)" if k == "4" else ""))
        ks = list(range(1, len(r["f" + k])))
        ax[1].semilogy(ks, [ab(c) for c in r["f" + k][1:]], "o-", ms=2.5, color=farbe, label=f"Rang {k} (M={M})")
    ax[0].set_xlabel("θ"); ax[0].set_title("Welle Φ_r(θ) − μ_r, Basen (e, 2)"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("Mode k"); ax[1].set_title("|Fourier-Koeffizient a_k|"); ax[1].legend(fontsize=8); ax[1].set_xlim(0, 65)
    fig.tight_layout(); fig.savefig(H / "wellen_raenge.png", dpi=130)
    print("Bild: wellen_raenge.png")
except Exception as e:
    print("kein Bild:", e)
