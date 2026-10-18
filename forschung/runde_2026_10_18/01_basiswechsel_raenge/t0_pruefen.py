"""Grundprüfungen der Bausteine (lib.py) gegen bekannte Werte und den Motor."""
import time
import mpmath as mp
import lib

mp.mp.dps = 100
t = time.time()
Te = lib.Tet("exp(1)")
T2 = lib.Tet("2")
print(f"Aufbau {time.time()-t:.1f}s")
for T in (Te, T2):
    print(T.name, "T(0)-1", mp.nstr(T.T(0) - 1, 3), "T(1)-b", mp.nstr(T.T(1) - T.b, 3), "T(-1)", mp.nstr(T.T(-1), 3))
    print("  L", mp.nstr(T.L, 30), "lam", mp.nstr(T.lam, 30), " g0", mp.nstr(T.g[0], 3))
    print("  slog(T(0.37))-0.37", mp.nstr(T.slog(T.T(mp.mpf("0.37"))) - mp.mpf("0.37"), 3))
    print("  P(0)-1", mp.nstr(T.P(0) - 1, 3), " P(1)-b", mp.nstr(T.P(1) - T.b, 3), " P(-1)", mp.nstr(T.P(-1), 3), " P(2)-T(b)", mp.nstr(T.P(2) - T.T(T.b), 3))
    print("  P(0.5)", mp.nstr(T.P(mp.mpf(0.5)), 30), " P(1.5)", mp.nstr(T.P(mp.mpf(1.5)), 30))
    print("  spog(P(0.3))-0.3", mp.nstr(T.spog(T.P(mp.mpf("0.3"))) - mp.mpf("0.3"), 3))
    print("  P(1.3) - T(P(0.3))", mp.nstr(T.P(mp.mpf("1.3")) - T.T(T.P(mp.mpf("0.3"))), 3))
# Rang 4: mu_{e,2}
t = time.time()
M = 32
v = [lib.phi4(Te, T2, mp.mpf(j) / M) for j in range(M)]
print("Phi4(e,2) Mittel (M=32)", mp.nstr(mp.fsum(v) / M, 30), f" {time.time()-t:.1f}s")
print("  Konvergenz in N:", [mp.nstr(lib.phi4(Te, T2, mp.mpf("0.3"), N) - lib.phi4(Te, T2, mp.mpf("0.3"), 30), 3) for N in (2, 3, 4, 5, 8)])
print("  max-min", mp.nstr(max(v) - min(v), 10))
