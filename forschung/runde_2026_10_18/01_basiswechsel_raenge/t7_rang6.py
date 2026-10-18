"""Rang 6: Prüfungen und erste Werte von Phi6."""
import sys, time
import mpmath as mp
import lib, lib6
mp.mp.dps = 120
bb, dd = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("exp(1)", "2")
Tb, Td = lib.Tet(bb), lib.Tet(dd)
Hb, Hd = lib6.Hex(Tb), lib6.Hex(Td)
for Hx in (Hb, Hd):
    print(Hx.T.name, "hex(0)-1", mp.nstr(Hx.hex(0) - 1, 3), "hex(1)-b", mp.nstr(Hx.hex(1) - Hx.T.b, 3),
          "hex(-1)", mp.nstr(Hx.hex(-1), 3), "ihex(hex(0.3))-0.3", mp.nstr(Hx.ihex(Hx.hex(mp.mpf("0.3"))) - mp.mpf("0.3"), 3),
          "hex(0.5)", mp.nstr(Hx.hex(mp.mpf(0.5)), 20))
th = mp.mpf("0.3")
print("G5x Kontrolle:", mp.nstr(lib6.G5x(Tb, Td, 5 + th) - 5 - th - lib.phi5(Tb, Td, th), 3),
      mp.nstr(lib6.G5x(Tb, Td, 1 + th) - lib.G5_endlich(Tb, Td, th, 1) - 1 - th, 3))
for n in (0, 1):
    print("endlich n", n, mp.nstr(lib6.G6_endlich(Tb, Td, Hb, Hd, th, n), 25))
for X in (30, 40):
    t = time.time()
    print("Phi6(0.3) X", X, mp.nstr(lib6.phi6(Tb, Td, Hb, Hd, th, X), 30), f"{time.time()-t:.1f}s")
print("Periodizität:", mp.nstr(lib6.phi6(Tb, Td, Hb, Hd, mp.mpf(0)) - lib6.phi6(Tb, Td, Hb, Hd, mp.mpf(1)), 3))
