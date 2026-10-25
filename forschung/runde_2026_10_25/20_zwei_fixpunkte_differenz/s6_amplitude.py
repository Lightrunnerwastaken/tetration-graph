"""s6: Amplitude C_p(b) = |c_1(p)| e^{S_a} und Phase von c_1 über weite Basisbereiche, mit Gegenprobe
(dps, N) gegen (dps+40, N+8). Aufruf: python s6_amplitude.py N K t1:dps1 t2:dps2 ...
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

N, K = int(sys.argv[1]), int(sys.argv[2])
print("t, b, lam_a, S_a, C_p (Einst.1), C_p (Einst.2), rel.Abw., Phase(c1 e^{2πi x1})/π")
for arg in sys.argv[3:]:
    t, dps = arg.split(":")
    out = []
    for d, n in [(int(dps), N), (int(dps) + 40, N + 8)]:
        mp.mp.dps = d
        P = Paar(mp.mpf(t), K=K)
        Sa = P.S_a()
        x1 = mp.re(P.alpha_a(mp.mpf(1)))
        vals = [P.p(mp.mpf(j) / n) for j in range(n)]
        c1 = mp.fsum(v * mp.expjpi(-2 * mp.mpf(j) / n) for j, v in enumerate(vals)) / n
        out.append((P, Sa, abs(c1) * mp.exp(Sa), mp.arg(c1 * mp.expjpi(2 * x1)) / mp.pi))
    mp.mp.dps = 20
    (P, Sa, C1, ph), (_, _, C2, _) = out
    print(f"{t} {mp.nstr(P.b, 15)} {mp.nstr(P.lam_a, 10)} {mp.nstr(Sa, 10)} {mp.nstr(C1, 14)} {mp.nstr(C2, 14)} {mp.nstr(abs(C1/C2-1), 3)} {mp.nstr(ph, 10)}", flush=True)
