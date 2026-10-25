"""s1: Kontrolle √2 (|c1| = 3.626e-25 bekannt) für p und q, dann Modenfolge bei kleinem t."""
import sys, time
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

t0 = time.time()
mp.mp.dps = 60
P = Paar(1, K=40)
for name, fun in [("p", P.p), ("q", P.q)]:
    cs = P.fourier(fun, 32)
    print(name, "c0", mp.nstr(cs[0], 12), "|c1|", mp.nstr(abs(cs[1]), 12), "|c1|e^S_a", mp.nstr(abs(cs[1]) * mp.exp(P.S_a()), 10))
print("zeit", time.time() - t0)
