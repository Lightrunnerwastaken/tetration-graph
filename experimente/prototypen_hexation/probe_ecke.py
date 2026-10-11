"""Höhe der Ecken der reellen Kurve in den W-Koordinaten des Fork-Verfahrens (W = α_p, Abel-Funktion am Fixpunkt p).
 θ(W) = Σ c_n e^{2πisnW} konvergiert nur oberhalb der Ecken; das Fork-Verfahren nimmt als Abklinghöhe 0.6·(htop − hmin) an.
 Ecke: W-Bild des Wertes L_T (bzw. P(L_T)), von oben angelaufen; Weg w = w0 + iy, y geometrisch von 0.8 bis 1e-12, knwalk macht die Zweige stetig.
 Aufruf: probe_ecke.py basis dps L_T
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, dps, LT = sys.argv[1], int(sys.argv[2]), sys.argv[3]
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
WEG = ("if(1, kys=vector(120,j,0.8*(1e-12/0.8)^((j-1)/119)); kw0={w0}; kpe=vector(120,j,kw0+I*kys[j]); kWe=knwalk(kpe); "
       "khe=vector(120,j,kns*imag(kWe[j])); khe[120])")
try:
    kopf = gp.eval_batch(basis, ["hexinit(sexp(1),0,1)", "hexerr", "knc", "knr", "kns", "knpar[1]", "knpar[2]", "knpar[3]", "knpar[4]", "knpar[5]",
                                 f"pent({LT})", "knres"])
    namen = ["Stellen", "knc", "knr", "kns", "K", "N", "Bogen/Im p", "hmin", "htop", "P(L_T)", "knres"]
    for n, v in zip(namen, kopf[1:]):
        print(f"  {n:11s} = {mp.nstr(v, 15)}", flush=True)
    hmin, htop, N = kopf[8].real, kopf[9].real, kopf[6].real
    hcut = hmin + mp.mpf("0.6") * (htop - hmin)
    for name, w0 in (("L_T", LT), ("P(L_T)", mp.nstr(kopf[10].real, dps))):
        r = gp.eval_batch(basis, [WEG.format(w0=w0), "real(kWe[120])", "imag(kWe[120])", "khe[60]", "khe[100]", "khe[119]"])
        hs = r[0].real
        print(f"  Ecke bei w = {name}: W = {mp.nstr(mp.mpc(r[1].real, r[2].real), 12)}, s·Im W = {mp.nstr(hs, 10)} (Weg: {mp.nstr(r[3].real, 6)}, "
              f"{mp.nstr(r[4].real, 6)}, {mp.nstr(r[5].real, 8)})", flush=True)
    rand_ann, rand_echt = hcut - hmin, hcut - hs
    print(f"  hcut = {mp.nstr(hcut, 8)}; angenommene Abklinghöhe hcut − hmin = {mp.nstr(rand_ann, 6)}, echte hcut − h_Ecke = {mp.nstr(rand_echt, 6)}")
    print(f"  Abbruchfehler der θ-Reihe am tiefsten θ-Punkt ~ e^(−2π N (hcut − h_Ecke)) = 1e{mp.nstr(-2 * mp.pi * N * rand_echt / mp.log(10), 4)}"
          f"  (angenommen 1e{mp.nstr(-2 * mp.pi * N * rand_ann / mp.log(10), 4)})")
finally:
    gp.close()
print(f"{time.time() - t0:.1f} s")
