"""Rang 10 (Basis 2): Kneser-Dekation des Motors, Fixpunkt kp10, x10* = H_10^{-1}(p_8), A10 = H_10'(x10*).
Aufruf: python p3_rang10.py DPS"""
import sys, json, time
from pathlib import Path
import mpmath as mp
import motor
DPS = int(sys.argv[1]); B = "2"; mp.mp.dps = DPS + 10
t0 = time.time()
r = motor.roh(B, ["precision(hyperdigits(9),10)", "hyp[9]", "hylam[9]", "hyc[9]"], DPS)
p8 = mp.re(mp.mpc(r[1]))
print("Rang 9:", r, f"[{time.time()-t0:.0f}s]", flush=True)
r = motor.roh(B, ["precision(hyperdigits(10),10)", "kp", "klam", "hyper(10,0.5)"], DPS)
print("Rang 10:", r, f"[{time.time()-t0:.0f}s]", flush=True)
x10 = mp.mpc(motor.roh(B, [f"ihyper(10,{motor.gpzahl(p8, DPS + 5)})"], DPS)[0])
h = mp.mpf(10) ** (-5)
w = [mp.mpc(v) for v in motor.roh(B, [f"hyper(10,{motor.gpzahl(x10 + h, DPS + 5)})", f"hyper(10,{motor.gpzahl(x10 - h, DPS + 5)})", f"hyper(10,{motor.gpzahl(x10, DPS + 5)})"], DPS)]
A10 = (w[0] - w[1]) / (2 * h)
print(f"x10* = {mp.nstr(x10, 18)}  H10(x10*) - p8 = {mp.nstr(w[2] - p8, 5)}  A10 = {mp.nstr(A10, 12)}  [{time.time()-t0:.0f}s]", flush=True)
json.dump({"p8": str(p8), "kp10": str(r[1]), "klam10": str(r[2]), "digits10": str(r[0]), "x10": str(x10), "A10": str(A10)},
          open(Path(__file__).resolve().parent / "daten" / f"rang10_{DPS}.json", "w"), indent=1)
motor.schliessen()
