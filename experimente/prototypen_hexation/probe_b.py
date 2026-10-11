import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 40
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=int(sys.argv[2]))
t0 = time.time()
try:
    ex = ["hexinit(sexp(1),0,1)", "hexreg", "hexerr", "kp", "klam", "knc", "kntl", "kntr"] + [f"hex({x})" for x in sys.argv[3].split(",")]
    for e, w in zip(ex, gp.eval_batch(sys.argv[1], ex)):
        print(f"{e:22s} {mp.nstr(w, 18)}", flush=True)
finally:
    gp.close()
print(f"{time.time() - t0:.1f} s")
