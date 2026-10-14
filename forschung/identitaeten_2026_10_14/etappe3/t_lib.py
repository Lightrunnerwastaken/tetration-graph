import mpmath as mp, time
mp.mp.dps = 60
from lib3 import Basis, ThetaReihe, lade_theta
from k3 import sexp
B = Basis(60)
print("L", B.L)
s = mp.mpc(0.3, 0.7)
w = B.F(s); print("h(F(s))-s", mp.nstr(B.h(w) - s, 5), "  (mod 2pi i/L =", mp.nstr(2j*mp.pi/B.L,8), ")")
print("Psi check", mp.nstr(B.Psi(B.L*mp.mpc(0.5,0.2)) - mp.exp(B.Psi(mp.mpc(0.5,0.2))), 5))
T = ThetaReihe(lade_theta(), B)
for z in [mp.mpc(0.3,0.2), mp.mpc(-0.4,0.15), mp.mpc(0.1,0.05)]:
    t=time.time(); a = T.sexp(z); tt=time.time()-t
    b = sexp([z], 80)[0]
    print(z, mp.nstr(a, 30), " diff", mp.nstr(abs(a-b), 3), f"{tt:.2f}s")
