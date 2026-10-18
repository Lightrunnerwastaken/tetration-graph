"""Gegenproben Rang 5: Periodizität von Phi5, Motor-Pentation (pentx/ipentx) gegen lib, endliche n."""
import mpmath as mp
import lib
import motor

mp.mp.dps = 150
Te, T2 = lib.Tet("exp(1)"), lib.Tet("2")
# 1) Periodizität und Stetigkeit
for th in (mp.mpf(0), mp.mpf(1), mp.mpf(2)):
    print("Phi5(e,2)", th, mp.nstr(lib.phi5(Te, T2, th), 40))
e = mp.mpf(10) ** -12
print("Sprung bei 1:", mp.nstr(lib.phi5(Te, T2, 1 - e) - lib.phi5(Te, T2, 1 + e), 5))
# 2) endliche n: G(n+th) - n - th
th = mp.mpf("0.3")
for n in range(0, 3):
    print("n", n, mp.nstr(lib.G5_endlich(Te, T2, th, n), 30))
print("Limes", mp.nstr(lib.phi5(Te, T2, th), 30))
# 3) Motor: pentx/ipentx bei dps 60
p2 = motor.werte("2", ["pentx(0.3)", "pentx(1.3)", "pentx(2.3)"], 60)
pe = motor.werte("exp(1)", [f"ipentx({motor.gpzahl(mp.re(p2[1]), 70)})", f"ipentx({motor.gpzahl(mp.re(p2[2]), 70)})"], 60)
print("P_2(0.3) Motor - lib", mp.nstr(mp.re(p2[0]) - T2.P(th), 3))
print("P_2(2.3) Motor - lib", mp.nstr(mp.re(p2[2]) - T2.P(2 + th), 3))
print("G(1.3)-1.3 Motor", mp.nstr(mp.re(pe[0]) - mp.mpf("1.3"), 30), " lib", mp.nstr(lib.G5_endlich(Te, T2, th, 1), 30))
print("G(2.3)-2.3 Motor", mp.nstr(mp.re(pe[1]) - mp.mpf("2.3"), 30), " lib", mp.nstr(lib.G5_endlich(Te, T2, th, 2), 30))
motor.schliessen()
