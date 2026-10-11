"""Strang C, c06: Wo konvergiert p_N ausserhalb von D? Cauchy-Test |p_N(x) - p_194(x)| an Testpunkten, die zwei
Literatur-Hypothesen trennen:
  (H1) Kreis um a0 = p(0) mit Radius |L - a0|   (p_N wie Taylor-Polynom von f_K an der Stelle a0)
  (H2) Walsh-Lemniskate |x (x - a0)| < |L (L - a0)|   (wie Hermite-Interpolation an den Knoten 0 und a0)
  (D)  |x| < |L|
Nur p_N aus gemeinsam/polynome_bis_194.json, kein Kneser-Motor.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, fabs, polyval

mp.dps = 60
HIER = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HIER, "..", "gemeinsam", "polynome_bis_194.json"), encoding="utf-8"))
L = mpc("0.31813150520476413", "1.3372357014306895")
a0 = mpf(P["194"][0])
rho2 = abs(L) * abs(L - a0)


def wert(N, x):
    c = [mpf(t) for t in P[str(N)]]
    return polyval(c[::-1], x)


punkte = [mpf("1.6"), mpf("1.8"), mpf("2.0"), mpc("0.5", "1.3"), mpc("0.5", "-1.3"), mpc("1.4", "0.9"),
          mpc("-0.6", "1.2"), mpf("-1.3"), mpc("1.0", "1.1")]
ref = {x: wert(194, x) for x in punkte}
print(f"a0 = {float(a0):.10f}, |L| = {float(abs(L)):.6f}, |L-a0| = {float(abs(L-a0)):.6f}, rho^2 = {float(rho2):.6f}")
for x in punkte:
    h1 = abs(x - a0) < abs(L - a0)
    h2 = abs(x) * abs(x - a0) < rho2
    hd = abs(x) < abs(L)
    zeile = []
    for N in (40, 80, 120, 160, 180):
        zeile.append(f"{float(abs(wert(N, x) - ref[x])):.1e}")
    print(f"x={complex(x):.2f}: |x-a0|={float(abs(x-a0)):.3f} H1={'ja' if h1 else 'nein'} "
          f"|x(x-a0)|={float(abs(x)*abs(x-a0)):.3f} H2={'ja' if h2 else 'nein'} D={'ja' if hd else 'nein'} | "
          f"|p_N - p_194| N=40,80,120,160,180: {' '.join(zeile)}")
