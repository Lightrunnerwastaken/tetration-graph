"""b26: Ist der Fehlerkern der J-Projektion separabel?  R_N(z,ζ) := 1 − (ζ − z)·Π_N k_ζ(z) (Polynom in z, R_N(ζ,ζ) = 1).
Bei Hermite-Interpolation wäre R = ω(z)/ω(ζ). Test: Kreuzverhältnis R(z,ζ)R(z',ζ')/(R(z,ζ')R(z',ζ)) ≈ 1 ?"""
import sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
N = 60
pw, g, Ji = baue(N, J=N)
zetas = [L * mp.mpf("1.1"), mp.mpc(-1.5, 0.3), mp.mpf("2.0"), mp.mpc(0.3, -1.6)]
zs = [mp.mpf("0.3"), mp.mpc(-0.8, 0.5), mp.mpc(1.0, 0.7), L * mp.mpf("0.98")]
Pz = {}
for zeta in zetas:
    P = loese(Ji, rhs_pol(N, pw, g, zeta))
    Pz[zeta] = P
R = lambda z, zeta: 1 - (zeta - z) * horner(Pz[zeta], z)
for i in range(len(zs) - 1):
    for j in range(len(zetas) - 1):
        z, z2, ze, ze2 = zs[i], zs[i + 1], zetas[j], zetas[j + 1]
        cr = R(z, ze) * R(z2, ze2) / (R(z, ze2) * R(z2, ze))
        print(f"z={mp.nstr(z,3)}, z'={mp.nstr(z2,3)}, ζ={mp.nstr(ze,3)}, ζ'={mp.nstr(ze2,3)}: Kreuzverhältnis = {mp.nstr(cr, 6)}")
