"""Strang A: Präzisionsprüfung Π_N f_K bei N = 800, 1000: 0.72N+200 Bit (gespeichert) gegen N+250 Bit (Vorschlag Leitung: N+200)."""
import sys, time
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from lin_kern import kneser_koeff, lade, pi_N
mp.mp.dps = 60
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
P = lade()
a = kneser_koeff(1001)
for N in [int(x) for x in sys.argv[1].split(",")]:
    t0 = time.time()
    b2 = pi_N(a, N, N + 250)
    b1 = [mp.mpf(x) for x in P[N]["b"]]
    b2 = [mp.mpf(x.mid().str(60, radius=False)) for x in b2]
    d = max(abs(x - y) * RL ** k for k, (x, y) in enumerate(zip(b1, b2)))
    hv = lambda c, z: mp.fsum(x * z ** k for k, x in enumerate(c))
    print(f"N = {N}: max_k |Δb_k||L|^k = {mp.nstr(d, 3)}; |ΔP(L)| = {mp.nstr(abs(hv(b1, L) - hv(b2, L)), 3)} "
          f"(|P(L)−L| = {mp.nstr(abs(hv(b1, L) - L), 4)}); {time.time() - t0:.0f} s", flush=True)
