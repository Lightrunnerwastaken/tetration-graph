"""Rang 6: Differenzenquotienten von Phi6(e,2) über viele Skalen; Ebenen V_j = P_d^j(hex_d(th0))."""
import mpmath as mp
import lib, lib6
mp.mp.dps = 120
Tb, Td = lib.Tet("exp(1)"), lib.Tet("2")
Hb, Hd = lib6.Hex(Tb), lib6.Hex(Td)
for th0 in ["0.05", "0.3"]:
    th0 = mp.mpf(th0)
    V = [Hd.hex(th0)]
    while V[-1] < 3:
        V.append(Td.P(V[-1]))
    print("th0", th0, "V_j =", [mp.nstr(v, 6) for v in V])
    prev = None
    for k in list(range(2, 41, 2)) + [50, 60]:
        h = mp.mpf(2) ** -k
        D = (lib6.phi6(Tb, Td, Hb, Hd, th0 + h) - lib6.phi6(Tb, Td, Hb, Hd, th0 - h)) / (2 * h)
        print(f"  h=2^-{k:<3d} D = {mp.nstr(D, 15)}" + (f"   Aenderung {mp.nstr(D - prev, 4)}" if prev is not None else ""), flush=True)
        prev = D
