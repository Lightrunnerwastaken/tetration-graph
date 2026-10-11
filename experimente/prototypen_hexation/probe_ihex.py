import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 40
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", "real(kp)", "imag(kp)", "knc", "kngood"])
    req, imq = k[1].real, k[2].real
    print("knc", mp.nstr(k[3], 8), "kngood", mp.nstr(k[4], 8), flush=True)
    for f in ("0", "0.3", "0.6", "0.8", "0.9", "0.95"):
        t0 = time.time()
        z = f"({mp.nstr(req, 30)}+I*{f}*{mp.nstr(imq, 30)})"
        w = gp.eval_batch("2", [f"iferr(ihex({z}),E,10^9)", f"iferr(hex(ihex({z}))-{z},E,10^9)"])
        print(f"t = {f} Im q: ihex = {mp.nstr(w[0], 12)}, Probe {mp.nstr(abs(w[1]), 3)}, {time.time() - t0:.1f} s", flush=True)
finally:
    gp.close()
