import sys, numpy as np
from flint import acb
N = int(sys.argv[1])
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", f"N = {N}"))
# y(t) = Σ_k y_k t^{-k} unskaliert: Z = D y (J^T Z = D u) ⇒ ⟨Z, g⟩ Paarung; y(x) als Funktional-Laurent: Y(x) = Σ_j Z_j x^{-j-1}
def Ylaurent(x):
    w = 1/acb(x.real, x.imag); s = acb(0)
    for v in Zl[::-1]: s = s*w + v
    s = s*w
    m = abs(s).mid(); return float(m.log())/N
Fp = np.polynomial.Polynomial(coef)
for xr in (1.05, 1.1, 1.2, 1.25, 1.3, 1.33):
    t = complex(xr*rho, 1e-9)
    r = (Fp - t).roots(); h = r[np.argmin(abs(r))]
    print(f"t = {xr:.2f}ρ: Φ_T = {T(t):+.4f};  h(t) = {abs(h)/rho:.4f}ρ;  Φ_y(h(t)) = {Ylaurent(h):+.4f};  Φ_y(t) = {Ylaurent(t):+.4f}", flush=True)
