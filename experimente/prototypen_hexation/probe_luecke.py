"""Wo unterscheiden sich die beiden Lösungen (zwei Bögen) des Fork-Verfahrens? Verdacht: in den Lücken des Bogens nahe der reellen Achse
 (knarc deckt nur die Winkel 0.05 … π − 0.05 ab). Ruft knplan/knsolve nur auf (die Fork-Datei bleibt unverändert).
 Aufruf: probe_luecke.py basis dps
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, dps = sys.argv[1], int(sys.argv[2])
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
LOES = ("if(1, kr1=knpar[3]*knr; kr2=knpar[10]*knr; kK=knpar[1]; kN=knpar[2]; kp1=knplan(kr1, precis-3); kp2=knplan(kr2, precis-3); "
        "kS1=knsolve(kr1, kK, kN, kp1[5]+0.6*(kp1[6]-kp1[5]), kp1[5]+0.85*(kp1[6]-kp1[5])); "
        "kS2=knsolve(kr2, kK, kN, kp2[5]+0.6*(kp2[6]-kp2[5]), kp2[5]+0.85*(kp2[6]-kp2[5])); "
        "kd=(u->abs(subst(kS1[1],x,(u-knc)/kS1[4])-subst(kS2[1],x,(u-knc)/kS2[4]))); kS1[5])")
try:
    kopf = gp.eval_batch(basis, ["hexinit(sexp(1),0,1)", "knc", "knr", "knpar[1]", "knpar[2]", "knpar[3]", "knpar[10]", LOES, "kS2[5]"])
    print(f"Basis {basis}, dps {dps}: {time.time() - t0:.1f} s; Stellen {mp.nstr(kopf[0], 5)}, c = {mp.nstr(kopf[1], 8)}, Im p = {mp.nstr(kopf[2], 8)}, "
          f"K = {int(kopf[3].real)}, N = {int(kopf[4].real)}, Bögen {mp.nstr(kopf[5], 3)} / {mp.nstr(kopf[6], 3)} Im p; Reste {mp.nstr(kopf[7], 3)}, {mp.nstr(kopf[8], 3)}", flush=True)
    knr, r1 = kopf[2].real, kopf[5].real * kopf[2].real
    winkel = [0, 0.005, 0.01, 0.02, 0.035, 0.05, 0.08, 0.15, 0.4, 0.8, 1.2, 1.5708]
    zeilen = []
    for f in (0.5, 0.9, 1.0):
        ex = []
        for phi in winkel:
            for seite in (1, -1):
                ex.append(f"kd(knc + {f}*kr1*exp(I*({phi} if {seite} > 0 else Pi - {phi})))".replace(f"({phi} if {seite} > 0 else Pi - {phi})",
                                                                                                 f"{phi}" if seite > 0 else f"Pi-{phi}"))
        werte = [v.real for v in gp.eval_batch(basis, ex)]
        zeilen.append((f, werte))
finally:
    gp.close()
print("Unterschied |A₁ − A₂| auf Kreisen um c (Radius f · Bogen₁), Winkel φ ab der rechten (φ) bzw. linken (π − φ) Achse:")
print("   φ        " + "   ".join(f"f={f}: rechts      links " for f, _ in zeilen))
for i, phi in enumerate(winkel):
    print(f"  {phi:7.4f}   " + "   ".join(f"{mp.nstr(w[2 * i], 3):>10s} {mp.nstr(w[2 * i + 1], 3):>10s}" for _, w in zeilen))
