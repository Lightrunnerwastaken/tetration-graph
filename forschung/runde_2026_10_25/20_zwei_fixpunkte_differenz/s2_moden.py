"""s2: Modenfolge c_k von p und q bei kleinem t (langsamer Abfall), Betrag und Phase.
r_k = c_k e^{S_a k}; Phase gegen -2πk x1 mit x1 = Re α_a(1) (Kandidat: Singularität von p bei α_a(1) = x1 - iτ_a).
Aufruf: python s2_moden.py t dps N K
"""
import sys, time
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

t, dps, N, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
mp.mp.dps = dps
t0 = time.time()
P = Paar(mp.mpf(t), K=K)
Sa, SA = P.S_a(), P.S_A()
al1 = P.alpha_a(mp.mpf(1))
x1 = mp.re(al1)
print(f"t={t} b={mp.nstr(P.b, 15)} a={mp.nstr(P.a, 12)} A={mp.nstr(P.A, 12)} lam_a={mp.nstr(P.lam_a, 10)} lam_A={mp.nstr(P.lam_A, 10)}")
print(f"S_a={mp.nstr(Sa, 12)} S_A={mp.nstr(SA, 12)}  alpha_a(1)={mp.nstr(al1, 15)}  tau_a={mp.nstr(mp.pi/abs(P.Lla), 12)}")
for name, fun in [("p", P.p), ("q", P.q)]:
    cs = P.fourier(fun, N)
    print(f"--- {name}: k, |c_k|, |r_k|=|c_k|e^(S_a k), |r_k|*k, arg(c_k e^(2πik x1))/π")
    for k in range(1, N // 2):
        r = cs[k] * mp.exp(Sa * k)
        ph = mp.arg(cs[k] * mp.expjpi(2 * k * x1)) / mp.pi
        print(f"{k:3d} {mp.nstr(abs(cs[k]), 8):>14} {mp.nstr(abs(r), 12):>18} {mp.nstr(abs(r)*k, 12):>18} {mp.nstr(ph, 10):>14}")
    sys.stdout.flush()
print("zeit", time.time() - t0)
