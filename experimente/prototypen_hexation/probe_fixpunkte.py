import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe
mp.mp.dps = 50
basis, start = sys.argv[1], sys.argv[2]
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
NEWTON = ("if(1, kz={s}; for(i=1,60, kf=pent(kz)-kz; kd=(pent(kz+1e-12)-pent(kz-1e-12))/2e-12-1; kz=kz-kf/kd; if(abs(kf)<1e-25,break)); kz)")
try:
    w = gp.eval_batch(basis, ["hexinit(sexp(1),0,1)", NEWTON.format(s=start), "abs(pent(kz)-kz)", "(pent(kz+1e-12)-pent(kz-1e-12))/2e-12", "kp", "knc"])
finally:
    gp.close()
z, res, mu, kp, knc = w[1], w[2], w[3], w[4], w[5]
print(f"Basis {basis}: Fixpunkt {mp.nstr(z, 15)}, |P(z) − z| = {mp.nstr(res, 3)}, P'(z) = {mp.nstr(mu, 10)}, |P'| = {mp.nstr(abs(mu), 8)}")
print(f"  abstossender Fixpunkt q = {mp.nstr(kp, 12)}; Abstand von c = Re q: {mp.nstr(abs(z - mp.re(kp)), 6)}, von knc = {mp.nstr(knc.real, 6)}: {mp.nstr(abs(z - knc.real), 6)}")
