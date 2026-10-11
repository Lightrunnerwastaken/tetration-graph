import sys, numpy as np
from flint import acb
N = int(sys.argv[1])
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", f"N = {N}"))
def Ylaurent(x):
    w = 1/acb(x.real, x.imag); s = acb(0)
    for v in Zl[::-1]: s = s*w + v
    s = s*w; m = abs(s).mid(); return float(m.log())/N
Fp = np.polynomial.Polynomial(coef)
argL = np.angle(L)
# Halbmond-Richtung bei L: zwischen Γout-Tangente (argL − 90°) und A-Bogen (argL − 68.9°): nimm argL − 80°
for psi_rel in (-80.0, -85.0):
    print(f"Richtung argL{psi_rel:+.0f}° (im Halbmond nahe L):")
    for r in (0.003, 0.01, 0.03, 0.06, 0.1, 0.2, 0.3):
        t = L + r*np.exp(1j*(argL + np.deg2rad(psi_rel)))
        rr_ = (Fp - t).roots(); h = rr_[np.argmin(abs(rr_ - L))]
        print(f"   r={r}: Φ_T={T(t):+.4f}  Φ_y(h)={Ylaurent(h):+.4f}  Φ_y(t)={Ylaurent(t):+.4f}  |h|/ρ={abs(h)/rho:.4f}", flush=True)
