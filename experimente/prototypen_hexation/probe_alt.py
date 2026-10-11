import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 40
basis, dps = sys.argv[1], int(sys.argv[2])
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
try:
    w = gp.eval_batch(basis, ["hexinit(sexp(1),0,1)", "if(1, knc=real(kp); kpl=knplan(0.6*knr, precis-3); 0)", "kpl[2]", "kpl[3]", "kpl[5]", "kpl[6]", "kpl[7]", "knr"])
finally:
    gp.close()
print(f"{basis} dps {dps}, knc = Re p, Bogen 0.6 Im p: K {mp.nstr(w[2],4)}, N {mp.nstr(w[3],4)}, hmin {mp.nstr(w[4],7)}, htop {mp.nstr(w[5],7)}, Reichweite {mp.nstr(w[6],6)} = {mp.nstr(w[6].real/w[7].real,4)} Im p")
