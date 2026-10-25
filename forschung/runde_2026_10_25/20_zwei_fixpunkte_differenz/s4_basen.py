"""s4: Abhängigkeit von der Basis. Für jedes t:
  C_p = |c_1(p)| e^{S_a},  Phase von c_1(p) gegen -2π Re α_a(1) - π/2,
  s0 (log-Singularität von q, Newton von der Schätzung aus c_{k+1}/c_k),  Gesetz für c_k(q) am grössten sicheren k.
Aufruf: python s4_basen.py dps N K kq t1 t2 ...
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

dps, N, K, kq = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
mp.mp.dps = dps
for t in sys.argv[5:]:
    P = Paar(mp.mpf(t), K=K)
    Sa = P.S_a()
    x1 = mp.re(P.alpha_a(mp.mpf(1)))
    cp = P.fourier(P.p, N)
    cq = P.fourier(P.q, N)
    ph = mp.arg(cp[1] * mp.expjpi(2 * x1)) / mp.pi
    line = (f"t={t} b={mp.nstr(P.b, 16)} lam_a={mp.nstr(P.lam_a, 10)} lam_A={mp.nstr(P.lam_A, 10)} S_a={mp.nstr(Sa, 12)} "
            f"|c1|={mp.nstr(abs(cp[1]), 12)} C_p=|c1|e^Sa={mp.nstr(abs(cp[1]) * mp.exp(Sa), 12)} phase1/pi={mp.nstr(ph, 8)}")
    print(line, flush=True)
    # Schätzung s0 aus Verhältnis
    rat = cq[kq + 1] / cq[kq]
    s_est = mp.log(rat) / (-2j * mp.pi)
    s_est = mp.mpc(mp.re(s_est) % 1, -abs(mp.im(s_est)))

    def Phi_a(z):
        n = 0
        while abs(z - P.a) > P.eps:
            z = P.f(z)
            n += 1
            if n > 20000:
                raise RuntimeError
        return P.sigma(z, P.a, P.qa, P.dqa) * P.lam_a ** (-n)
    try:
        s0 = mp.findroot(lambda s: Phi_a(P.T_A(s)), s_est)
        laws = []
        for k in range(1, kq + 2):
            pred = -mp.exp(-2j * mp.pi * k * s0) / (k * P.Lla)
            laws.append(mp.nstr(abs(cq[k] / pred - 1), 3))
        print(f"   s0={mp.nstr(s0, 18)} (Schätzung {mp.nstr(s_est, 6)})  2π|Im s0|-S_a={mp.nstr(2*mp.pi*abs(mp.im(s0)) - Sa, 12)}  "
              f"T_A(s0)={mp.nstr(P.T_A(s0), 10)}  |c_k(q)/Gesetz-1| k=1..{kq+1}: {laws}")
    except Exception as e:
        print("   s0 nicht gefunden", e, mp.nstr(s_est, 8))
    sys.stdout.flush()
