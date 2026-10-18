"""Rang 6: c6(V) = H6(V) - V nähert sich Phi5({V}) (Rang-5-Phase) für grosse V; Gewicht der Ebene 2 in Phi6."""
import mpmath as mp
import lib, lib6
mp.mp.dps = 120
Tb, Td = lib.Tet("exp(1)"), lib.Tet("2")
Hb, Hd = lib6.Hex(Tb), lib6.Hex(Td)
for V in ["2.3", "3.3", "4.3", "6.6", "10.5", "40.25"]:
    V = mp.mpf(V)
    c6 = lib6.H6(Tb, Td, V) - V
    p5 = lib.phi5(Tb, Td, V - mp.floor(V))
    print(f"V={float(V):6.2f}  c6(V) = {mp.nstr(c6, 15):>20}  Phi5(frac V) = {mp.nstr(p5, 15):>20}  Differenz = {mp.nstr(c6 - p5, 3)}")
