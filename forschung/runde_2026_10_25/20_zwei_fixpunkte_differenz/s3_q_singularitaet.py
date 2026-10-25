"""s3: Logarithmische Singularität von q(s) = α_a(T_A(s)) - s: Nullstelle s0 von Φ_a(T_A(s)) (also T_A(s0) Urbild von a),
dann Vergleich c_k(q) mit dem Gesetz  c_k ≈ -exp(-2πi k s0)/(k ln λ_a)   (Singularität in der unteren Halbebene).
T_A ist ganz (nur exp), Φ_a nur Vorwärtsiteration: keine Zweigwahl nötig.
Aufruf: python s3_q_singularitaet.py t dps N K re_start im_start
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

t, dps, N, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
s_start = mp.mpc(sys.argv[5], sys.argv[6])
mp.mp.dps = dps
P = Paar(mp.mpf(t), K=K)


def Phi_a(z):
    n = 0
    while abs(z - P.a) > P.eps:
        z = P.f(z)
        n += 1
        if n > 100000:
            raise RuntimeError("keine Konvergenz")
    u = P.sigma(z, P.a, P.qa, P.dqa)
    return u * P.lam_a ** (-n)


s0 = mp.findroot(lambda s: Phi_a(P.T_A(s)), s_start)
z0 = P.T_A(s0)
print(f"t={t} b={mp.nstr(P.b, 15)}  s0={mp.nstr(s0, 25)}")
print(f"  T_A(s0)={mp.nstr(z0, 20)}   (a={mp.nstr(P.a, 20)}, 2π/ln b={mp.nstr(2*mp.pi/P.c, 15)})")
print(f"  f(T_A(s0))={mp.nstr(P.f(z0), 15)}  T_A(s0+1)={mp.nstr(P.T_A(s0+1), 15)}")
print(f"  2π|Im s0|={mp.nstr(2*mp.pi*abs(mp.im(s0)), 15)}  S_a={mp.nstr(P.S_a(), 15)}  π/ln λ_A={mp.nstr(mp.pi/P.LlA, 12)}  τ_a={mp.nstr(mp.pi/abs(P.Lla), 12)}")
cs = P.fourier(P.q, N)
print("  k   |c_k(q)|        Verhältnis c_k / (-e^{-2πiks0}/(k ln λ_a))")
for k in range(1, N // 2):
    pred = -mp.exp(-2j * mp.pi * k * s0) / (k * P.Lla)
    r = cs[k] / pred
    print(f"{k:3d} {mp.nstr(abs(cs[k]), 8):>14}  {mp.nstr(r, 14)}")
