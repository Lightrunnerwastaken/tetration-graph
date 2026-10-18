"""Kozykel der Rang-5-Phase: Phi_{b,d}(th) = Phi_{c,d}(th) + Phi_{b,c}(th + Phi_{c,d}(th)), und inverse Gleichung."""
import mpmath as mp
import lib
mp.mp.dps = 120
T2, Te, T3 = lib.Tet("2"), lib.Tet("exp(1)"), lib.Tet("3")
for th in ["0.1", "0.4", "0.8"]:
    th = mp.mpf(th)
    pbd = lib.phi5(T2, T3, th)
    pcd = lib.phi5(Te, T3, th)
    pbc = lib.phi5(T2, Te, th + pcd)
    inv = lib.phi5(T3, T2, th + pbd)
    print(f"th={float(th)}: Phi(2,3) = {mp.nstr(pbd, 25)}  Kozykel-Rest = {mp.nstr(pbd - pcd - pbc, 3)}  inverse: {mp.nstr(inv + pbd, 3)}")
