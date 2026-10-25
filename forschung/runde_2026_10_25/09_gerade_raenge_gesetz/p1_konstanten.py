"""Konstanten der Motor-Familie zur Basis 2: p_k, lambda_k (k gerade), Kneser-Fixpunkte kp_r, x_r* = H_r^{-1}(p_{r-2}),
vorhergesagte Exponenten alpha_r = ln lambda_{r-4} / ln lambda_{r-2}.
Aufruf: python p1_konstanten.py DPS RMAX   (RMAX gerade, 8 oder 10)"""
import sys, json, time
from pathlib import Path
import mpmath as mp
import motor
DPS, RMAX = int(sys.argv[1]), int(sys.argv[2])
B = "2"
t0 = time.time()
out = {}
r = motor.roh(B, ["pentx(0.5)", "hexx(0.5)", "lambdat", "pcent", "psing", "precision(hyperdigits(7),10)"], DPS)
print("pent/hex ok", r, f"{time.time()-t0:.0f}s", flush=True)
out["lambda4"], out["pcent5"], out["psing6"] = str(r[2]), str(r[3]), str(r[4])
# p_4 = reeller Fixpunkt von sexp (Pentationsfixpunkt): aus ipent(-oo)?  der Motor nennt ihn in pentinit; hier per Newton
for R in range(8, RMAX + 1, 2):
    t = time.time()
    r = motor.roh(B, [f"precision(hyperdigits({R}),10)", "kp", "klam", f"hyp[{R-1}]", f"hylam[{R-1}]", f"hyc[{R-1}]"], DPS)
    print(f"Rang {R}: Stellen {r[0]}  kp={r[1]}  klam={r[2]}  p_{R-2}={r[3]}  lam_{R-2}={r[4]}  hyc={r[5]}  [{time.time()-t:.0f}s]", flush=True)
    p = r[3]
    xs = motor.roh(B, [f"ihyper({R},{motor.gpzahl(p, DPS + 5)})"], DPS)[0]
    print(f"   x_{R}* = ihyper({R}, p_{R-2}) = {xs}", flush=True)
    out[f"r{R}"] = {"digits": str(r[0]), "kp": str(r[1]), "klam": str(r[2]), "p": str(r[3]), "lam": str(r[4]), "hyc": str(r[5]), "xstar": str(xs)}
json.dump(out, open(Path(__file__).resolve().parent / "daten" / f"konstanten_{DPS}_{RMAX}.json", "w"), indent=1)
print(f"[{time.time()-t0:.0f}s] fertig", flush=True)
motor.schliessen()
