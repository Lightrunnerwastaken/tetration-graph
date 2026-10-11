"""Strang A, Aufgabe 3/4: Modengehalt von E_N = p_N − f_K in Abel-Koordinaten.

Liegt E_N im Kern der Linearisierung (E∘f + (f'∘f)E = 0), dann ist θ_N(w) := E_N(T(w)) / T'(w + 1/2) eine Funktion von w allein,
1-periodisch und ½-antiperiodisch: θ_N(w) = Σ_{n ungerade} c_n e^{2πinw}, und E_N = Σ c_n δ_n.
Gemessen auf waagrechten Strecken w = x + iy, x ∈ [x0, x0 + 1), y = 0, 0.25, 0.5 (z = T(w) innerhalb |z| < 1):
    C_n(y) = Mittel θ_N e^{−2πinx},  c_n = C_n(y)·e^{2πny}.
Prüfungen: gerade n (inkl. n = 0) ≈ 0; c_n unabhängig von y. Dann c_{−1}(N) (Mode singulär bei L) gegen N.
Modell für x = 1: p_N(1) − T(1/2) = −r_N(1) + 2 Re(c_{−1} δ_{−1}(1)) + …
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 75
ctx.prec = 300
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
H = mp.mpf(1) / 2

q = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
P = {int(k): [mp.mpf(x) for x in v] for k, v in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8")).items()}
# slog-Koeffizienten bis 1500 aus den Kreisdaten r = 1.34
kd = json.loads((DATEN / "kreis_r1.34_M4096.json").read_text(encoding="utf-8"))
M = kd["M"]
Ah = [mp.mpc(*x) for x in kd["A"]]
Av = Ah + [mp.conj(x) for x in Ah[1:-1][::-1]]
dd = acb.dft([acb(arb(mp.nstr(mp.re(x), 70)), arb(mp.nstr(mp.im(x), 70))) for x in Av])
r = mp.mpf("1.34")
Akoef = [mp.mpf(dd[k].real.mid().str(70, radius=False)) / M / r ** k for k in range(1500)]
Tq = mp.taylor  # noqa: F841


def reihe(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def abl(c, z):
    return reihe([k * c[k] for k in range(1, len(c))], z)


m = 64
X0 = mp.mpf("-1.6")
YS = [mp.mpf(0), mp.mpf("0.25"), mp.mpf("0.5")]
gp = FatouGP(dps=57, fatou_gp="fork")
s = gp.session("exp(1)")
linien = {}
for y in YS:
    ws = [mp.mpc(X0 + mp.mpf(j) / m, y) for j in range(m)]
    zs = s.sexp_batch(ws)
    if y == 0:
        zs = [mp.mpc(mp.re(z), 0) for z in zs]
    Tp = [abl(q, z) / abl(Akoef, z) for z in zs]          # T'(w + 1/2) = f'(z)/A'(z)
    fz = [reihe(q, z) for z in zs]
    linien[y] = (ws, zs, Tp, fz)
    print(f"Linie y = {y}: max|z| = {mp.nstr(max(abs(z) for z in zs), 4)}, "
          f"Prüfung T(w+½) = f(z): {mp.nstr(max(abs(a - b) for a, b in zip(s.sexp_batch([w + H for w in ws]), fz)), 3)}")
# δ_{−1}(1) = T'(A(1)+½)·e^{−2πiA(1)} = T'(1/2)·1  (A(1) = 0)
d1_1 = abl(q, mp.mpf(1)) / abl(Akoef, mp.mpf(1))
gp.close()
print(f"T'(1/2) = δ_(±1)(1) = {mp.nstr(d1_1, 20)}")

NS_MODE = [-5, -3, -2, -1, 0, 1, 2, 3, 5]
erg = {}
print("\n N   |c_-1|        arg c_-1   |c_0|/|c_-1|  |c_±2|/|c_-1|  |c_-3|/|c_-1|  y-Konsistenz(c_-1)  "
      "p_N(1)−T(½)     Modell −r_N(1)+2Re(c_-1 δ_-1(1))")
for N in list(range(10, 200, 5)) + [194]:
    if N > 194 or N in erg:
        continue
    pN = P[N]
    cs = {}
    for y in YS:
        ws, zs, Tp, fz = linien[y]
        th = [(reihe(pN, z) - f) / t for z, f, t in zip(zs, fz, Tp)]
        cs[y] = {n: mp.fsum(t * mp.expjpi(-2 * n * (mp.re(w) - X0)) for t, w in zip(th, ws)) / m
                 * mp.expjpi(-2 * n * X0) * mp.exp(2 * mp.pi * n * y) for n in NS_MODE}
    c = cs[YS[0]]
    kons = max(abs(cs[y][-1] - c[-1]) for y in YS) / abs(c[-1])
    wert = mp.fsum(pN) - mp.fsum(q)               # p_N(1) − f_K(1), f_K(1) = Σ a_k (Schwanz ab 1001 vernachlässigbar)
    rN1 = mp.fsum(q[N + 1:])
    modell = -rN1 + 2 * mp.re(c[-1] * d1_1)
    erg[N] = {"c": c, "kons": kons, "wert": wert, "modell": modell, "rN1": rN1}
    print(f"{N:3d}  {mp.nstr(abs(c[-1]), 5):>11}  {mp.nstr(mp.arg(c[-1]), 5):>8}  {mp.nstr(abs(c[0]) / abs(c[-1]), 2):>9}  "
          f"{mp.nstr(max(abs(c[2]), abs(c[-2])) / abs(c[-1]), 2):>9}  {mp.nstr(abs(c[-3]) / abs(c[-1]), 2):>9}  {mp.nstr(kons, 2):>9}  "
          f"{mp.nstr(wert, 5):>12}  {mp.nstr(modell, 5):>12}", flush=True)

print("\nlokale Exponenten von c_{-1}(N) ~ N^{−β}: Re β aus |c|, Im β aus arg (Soll 4πi/L = "
      f"{mp.nstr(4j * mp.pi / L, 6)}):")
Ns = sorted(erg)
for N1, N2 in zip(Ns[:-1], Ns[1:]):
    a, b = erg[N1]["c"][-1], erg[N2]["c"][-1]
    lq = mp.log(mp.mpf(N2) / N1)
    dphi = mp.arg(b / a)
    print(f"  {N1:3d}→{N2:3d}: Re β = {mp.nstr(-mp.log(abs(b) / abs(a)) / lq, 4)},  Im β = {mp.nstr(-dphi / lq, 4)}")
(DATEN / "abel_moden.json").write_text(json.dumps({str(N): {"c": {str(n): [mp.nstr(mp.re(v), 20), mp.nstr(mp.im(v), 20)]
                                                                   for n, v in e["c"].items()},
                                                             "kons": mp.nstr(e["kons"], 3), "wert": mp.nstr(e["wert"], 20),
                                                             "modell": mp.nstr(e["modell"], 20), "rN1": mp.nstr(e["rN1"], 20)}
                                                   for N, e in erg.items()}, indent=0), encoding="utf-8")
print("geschrieben", time.strftime("%H:%M"))
