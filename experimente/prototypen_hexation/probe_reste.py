"""Welche Gleichungen von knsolve tragen den Restfehler? θ-Gleichungen (oben, nahe q) gegen Funktionalgleichungen (unten, nahe der Achse).
 Rechnet die Bogenpunkte von knhexinit nach (gleiche knarc/knwalk-Aufrufe) und wertet die Lösung knA/kncc an ihnen aus.
 Aufruf: probe_reste.py basis dps
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, dps = sys.argv[1], int(sys.argv[2])
kraft = sys.argv[3] if len(sys.argv) > 3 else "0"
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
SETUP = ("if(1, kr=knpar[3]*knr; kpts=knarc(kr, ceil(0.8*(knpar[1]+1+2*(knpar[2]+1)))); kW=knwalk(kpts); "
         "kh=vector(#kW,j,kns*imag(kW[j])); khc=knpar[4]+0.6*(knpar[5]-knpar[4]); khc2=knpar[4]+0.85*(knpar[5]-knpar[4]); "
         "kRt=vector(#kpts); kRf=vector(#kpts); kG=vector(#kpts); "
         "for(j=1,#kpts, if(kh[j]>=khc, kRt[j]=abs(knser(kpts[j])-kW[j]-sum(n=1,#kncc,kncc[n]*knen(kW[j],n-1)))); "
         "if(kh[j]<khc2, kG[j]=if(real(kpts[j])>knc, ipent(kpts[j]), pent(kpts[j])); "
         "kRf[j]=abs(knser(kpts[j])-knser(kG[j])-if(real(kpts[j])>knc,1,-1)))); #kpts)")
try:
    w = gp.eval_batch(basis, [f"if(1, knforce={kraft}; knmaxk=2000; 0)", "hexinit(sexp(1),0,1)", "hexerr", "knres", SETUP, "knc", "knpar[3]*knr", "knpar[1]"])[1:]
    n = int(w[3].real)
    knc, rad, K = w[4].real, w[5].real, int(w[6].real)
    print(f"knc = {mp.nstr(knc, 8)}, Bogenradius {mp.nstr(rad, 6)}, K = {K}", flush=True)
    print(f"Basis {basis}, dps {dps}: {time.time() - t0:.1f} s, Stellen {mp.nstr(w[1], 5)}, knres {mp.nstr(w[2], 3)}, {n} Bogenpunkte", flush=True)
    ex = []
    for j in range(1, n + 1):
        ex += [f"real(kpts[{j}])", f"imag(kpts[{j}])", f"kh[{j}]", f"kRt[{j}]", f"kRf[{j}]", f"real(kG[{j}])", f"imag(kG[{j}])"]
    v = [x.real for x in gp.eval_batch(basis, ex)]
finally:
    gp.close()
zeilen = [v[7 * j:7 * j + 7] for j in range(n)]
c = None
print(" j   Winkel   s·Im W     Rest θ      Rest FE    |Bild − c|")
for j, (x, y, h, rt, rf, gx, gy) in enumerate(zeilen):
    if j % max(1, n // 40) == 0 or rt > 1e-27 or rf > 1e-27:
        ang = mp.atan2(y, x - zeilen[0][0] + (zeilen[0][0] - (zeilen[0][0] + zeilen[-1][0]) / 2) * 0) if False else None
        print(f"{j + 1:3d}  {mp.nstr(mp.mpc(x, y), 8):28s} h = {mp.nstr(h, 5):8s} θ: {mp.nstr(rt, 3):10s} FE: {mp.nstr(rf, 3):10s} Bild {mp.nstr(mp.mpc(gx, gy), 8)}")
mt = max(range(n), key=lambda j: zeilen[j][3])
mf = max(range(n), key=lambda j: zeilen[j][4])
fe = [(abs(mp.mpc(z[5], z[6]) - knc) / rad, z[4]) for z in zeilen if z[4] > 0]
fe.sort()
print("FE-Reste nach Lage des Bildpunkts (|g − c| / Bogenradius):")
for lo, hi in ((0, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, 1.2), (1.2, 2)):
    t = [r for q, r in fe if lo <= q < hi]
    if t:
        print(f"   {lo:.1f} … {hi:.1f}: {len(t):3d} Gleichungen, grösster Rest {mp.nstr(max(t), 3)}, Mittel {mp.nstr(sum(t) / len(t), 3)}")
print(f"grösster θ-Rest {mp.nstr(zeilen[mt][3], 3)} bei Punkt {mt + 1}, grösster FE-Rest {mp.nstr(zeilen[mf][4], 3)} bei Punkt {mf + 1}")
