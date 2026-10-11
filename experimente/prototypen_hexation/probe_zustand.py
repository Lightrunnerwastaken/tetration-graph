import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 60
basis, dps = sys.argv[1], int(sys.argv[2])
for lauf in range(2):
    gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
    t0 = time.time()
    try:
        ex = ["hexinit(sexp(1),0,1)", "hexerr", "knres", "kp", "knc", "knpar[1]", "knpar[2]", "knpar[3]", "default(realprecision)", "sexp(0.5)"]
        w1 = gp.eval_batch(basis, ex)
        w2 = gp.eval_batch(basis, ["kp", "knc", "hexerr", "default(realprecision)"])
    finally:
        gp.close()
    print(f"Lauf {lauf}: {time.time() - t0:.1f} s", flush=True)
    for e, v in zip(ex, w1):
        print(f"   {e:24s} {mp.nstr(v, 20)}")
    print("   zweiter Aufruf:", [mp.nstr(v, 12) for v in w2])
