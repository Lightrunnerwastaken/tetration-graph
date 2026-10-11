"""Kurze Probe: Kneser-Hexation aus dem gp-tetration-Fork (knhexinit, exp-079) — Zeit, Werte, θ-Konstante; Aufruf: probe_fork.py basis [dps]."""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis = sys.argv[1]
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
try:
    namen = ["hexinit(sexp(1),0,1)", "hexreg", "hexerr", "kp", "klam", "kns", "kncc[1]", "kncc[2]", "knpar[1]", "knpar[2]", "knpar[3]", "knpar[7]",
             "hex(0.5)", "hex(1.5)", "hex(1)-sexp(1)", "hex(-3)+2"]
    werte = gp.eval_batch(basis, namen)
    print(f"Basis {basis}, dps {dps}: {time.time() - t0:.1f} s", flush=True)
    for n, w in zip(namen, werte):
        print(f"  {n:22s} = {mp.nstr(w, 25)}", flush=True)
finally:
    gp.close()
