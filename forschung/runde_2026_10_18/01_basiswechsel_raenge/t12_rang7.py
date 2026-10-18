"""Rang 7: Prüfungen und erste Werte von Phi7."""
import sys, time
import mpmath as mp
import lib, lib6, lib7
mp.mp.dps = 120
bb, dd = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("exp(1)", "2")
Tb, Td = lib.Tet(bb), lib.Tet(dd)
Hb, Hd = lib6.Hex(Tb), lib6.Hex(Td)
Yb, Yd = lib7.Hyp7(Hb), lib7.Hyp7(Hd)
for Y in (Yb, Yd):
    print(Y.Hx.T.name, "Tschebyschow-Schwanz", mp.nstr(Y.schwanz, 3), "hyp(0)-1", mp.nstr(Y.hyp(0) - 1, 3), "hyp(1)-b", mp.nstr(Y.hyp(1) - Y.Hx.T.b, 3),
          "ihyp(hyp(0.3))-0.3", mp.nstr(Y.ihyp(Y.hyp(mp.mpf("0.3"))) - mp.mpf("0.3"), 3), "hyp(0.5)", mp.nstr(Y.hyp(mp.mpf(0.5)), 20))
th = mp.mpf("0.3")
for X in (20, 30):
    t = time.time()
    print("Phi7(0.3) X", X, mp.nstr(lib7.phi7(Tb, Td, Hb, Hd, Yb, Yd, th, X), 30), f"{time.time()-t:.1f}s", flush=True)
print("Periodizität:", mp.nstr(lib7.phi7(Tb, Td, Hb, Hd, Yb, Yd, mp.mpf(0)) - lib7.phi7(Tb, Td, Hb, Hd, Yb, Yd, mp.mpf(1)), 3))
