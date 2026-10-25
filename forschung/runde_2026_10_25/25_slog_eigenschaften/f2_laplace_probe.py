"""F2: Gegenprobe der exakten Laplace-Darstellung der Fourier-Koeffizienten der Welle P (Basis e, Kneser):
    a_n = -i * int_0^oo [ sexp(-iy) - sexp(-1-iy) - 1 ] e^{-2 pi n y} dy      (n >= 1)
Herleitung: Strecke [-1,0] nach unten ins Unendliche verschieben (sexp -> conj(L), e^{-2 pi i n u} faellt), Periodizitaet von e^{-2 pi i n u}.
sexp(-iy): Taylor um 0 (Radius 2) fuer y <= 1.2, sonst gp-Motor; sexp(-1-iy) = log(sexp(-iy)) (Hauptzweig, Im wird geprueft).
Zweite Einstellung: Gauss-Legendre mit 2 verschiedenen Knotenzahlen; Vergleich mit direkter Quadratur aus f1.
"""
import json, sys, time
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
from motor import werte, schliessen, gpzahl

mp.mp.dps = 30
C0 = [mp.mpf(s) for s in json.load(open(HIER / "daten" / "taylor_0_1_512_60.json"))][:200]
S_taylor = lambda z: mp.polyval(C0[::-1], z)

Y = 14  # e^{-2 pi Y} ~ 1e-38
def knoten(grad, teile):
    g = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(grad, mp.mp.prec)
    kn = []
    grenzen = [mp.mpf(0), mp.mpf("1.2")] + [mp.mpf("1.2") + (Y - mp.mpf("1.2")) * j / teile for j in range(1, teile + 1)]
    for lo, hi in zip(grenzen[:-1], grenzen[1:]):
        for x, w in g:
            kn.append(((lo + hi) / 2 + (hi - lo) / 2 * x, w * (hi - lo) / 2))
    return kn

def sexp_achse(ys):
    innen = [y for y in ys if y <= mp.mpf("1.2")]
    aussen = [y for y in ys if y > mp.mpf("1.2")]
    res = {}
    for y in innen:
        res[y] = S_taylor(mp.mpc(0, -y))
    if aussen:
        t = time.time()
        v = werte("exp(1)", [f"sexp({gpzahl(mp.mpc(0, -y), 40)})" for y in aussen], 30)
        print(f"  Motor: {len(aussen)} Werte in {time.time()-t:.1f}s", flush=True)
        for y, w in zip(aussen, v):
            res[y] = w
    return res

for grad, teile in [(4, 4), (5, 5)]:
    kn = knoten(grad, teile)
    print(f"Einstellung Gauss-Grad {grad}, {teile} Aussenteile: {len(kn)} Knoten", flush=True)
    sv = sexp_achse([y for y, w in kn])
    maxim = max(abs(mp.im(mp.log(sv[y]))) for y, w in kn)
    print("  max |Im log sexp(-iy)| =", mp.nstr(maxim, 8), "(< pi: Hauptzweig konsistent)")
    # Probe des Zweigs an der Nahtstelle: Taylor um -1 nicht noetig; Probe: |sexp(-iy) - conj(L)| bei y = Y
    L = mp.mpc(0.31813150520476413531, 1.3372357014306894089)
    L = mp.findroot(lambda z: mp.exp(z) - z, L)
    yM = max(y for y, w in kn)
    print("  |sexp(-iY) - conj(L)| bei Y =", mp.nstr(yM, 6), ":", mp.nstr(abs(sv[yM] - mp.conj(L)), 5))
    for n in [1, 2, 3, 5]:
        s = mp.fsum(w * (sv[y] - mp.log(sv[y]) - 1) * mp.exp(-2 * mp.pi * n * y) for y, w in kn)
        a = -1j * s
        print(f"  n={n}: a_n (Laplace) = {mp.nstr(a.real, 20)} + {mp.nstr(a.imag, 20)} i", flush=True)
schliessen()
