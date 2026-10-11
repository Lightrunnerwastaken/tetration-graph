import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 40
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
try:
    for ex in (["if(1, knforce=[0.6,0,0,0,0]; 0)", "hexinit(sexp(1),0,1)", "hexreg", "hexerr", "poldegree(knA)", "knc", "real(kp)"],
               ["if(1, knforce=0; 0)", "hexinit(sexp(1),0,1)", "hexreg", "hexerr", "poldegree(knA)", "knc"]):
        try:
            print([mp.nstr(v, 10) for v in gp.eval_batch("2", ex)], flush=True)
        except Exception as e:
            print("Fehler:", str(e)[:300], flush=True)
finally:
    gp.close()
