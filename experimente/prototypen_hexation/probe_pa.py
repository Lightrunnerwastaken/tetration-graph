import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 40
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=20)
t0 = time.time()
try:
    z = "(0.33+1.37*I)"
    w = gp.eval_batch("exp(1)", ["pentinit(sexp(1),1)", f"iferr(pent({z}),E,0)", f"iferr(pent(pent({z})),E,0)", f"iferr(pent(pent(pent({z}))),E,0)"])
    print([mp.nstr(v, 10) for v in w], f"{time.time() - t0:.1f} s", flush=True)
finally:
    gp.close()
