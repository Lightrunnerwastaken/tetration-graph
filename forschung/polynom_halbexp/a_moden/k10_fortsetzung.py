"""Strang A: p_N über N = 194 hinaus fortsetzen (Newton wie experimente/polynom_wurzel.py, Start p_{N−1} mit a_N = 0),
um die Vorhersage „Vorzeichenwechsel von p_N(1) − T(½) bei N ≈ 217“ (k9_konstante.py) blind zu prüfen.
Aufruf: python k10_fortsetzung.py NMAX [bits]. Ergebnis laufend in daten/polynome_ab_195.json.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Matura\tetration-graph\experimente")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import arb, ctx  # noqa: E402

from flint import arb_mat, arb_poly  # noqa: E402


def _schritt(a, N, tol=mp.mpf(10) ** -75):
    """Newton wie polynom_wurzel._schritt, aber Abbruch bei |d| < tol oder Stillstand (Rauschen ~ 1.63^N·2^−prec)."""
    alt = None
    for it in range(40):
        f = arb_poly(a)
        P = [arb_poly([1])]
        for _m in range(N):
            P.append((P[-1] * f).truncate(N + 1))
        g, fpf = arb_poly([0]), arb_poly([0])
        for m in range(N + 1):
            g += a[m] * P[m]
            if m:
                fpf += m * a[m] * P[m - 1]
        F, J, fak = arb_mat(N + 1, 1), arb_mat(N + 1, N + 1), arb(1)
        for k in range(N + 1):
            fak = fak * k if k else fak
            F[k, 0] = g[k] - 1 / fak
            for j in range(N + 1):
                J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
        d = J.solve(F)
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
        gr = max(abs(mp.mpf(d[i, 0].mid().str(30, radius=False))) for i in range(N + 1))
        if gr < tol or (alt is not None and gr > alt / 10 and gr < mp.mpf(10) ** -40):
            return a, gr, it
        alt = gr
    raise RuntimeError(f"Newton ohne Konvergenz bei N = {N}")

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
NMAX = int(sys.argv[1])
ctx.prec = int(sys.argv[2]) if len(sys.argv) > 2 else 600
mp.mp.dps = 80
T12 = mp.mpf("1.6463542337511945809719240315921145182")    # T_e(1/2) (Σ a_k, 40 Stellen)
pfad = DATEN / "polynome_ab_195.json"
erg = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {}
if erg:
    N0 = max(int(k) for k in erg)
    a = [arb(x) for x in erg[str(N0)]]
else:
    N0 = 194
    a = [arb(x) for x in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8"))["194"]]
t0 = time.time()
for N in range(N0 + 1, NMAX + 1):
    a, gr, it = _schritt(a + [arb(0)], N)
    w = mp.fsum(mp.mpf(c.mid().str(70, radius=False)) for c in a)
    erg[str(N)] = [c.mid().str(70, radius=False) for c in a]
    pfad.write_text(json.dumps(erg), encoding="utf-8")
    print(f"N = {N}: p_N(1) − T(½) = {mp.nstr(w - T12, 6)}   (Newton {it + 1} Schritte, letzte Korrektur {mp.nstr(gr, 2)}; {time.time() - t0:.0f} s)", flush=True)
