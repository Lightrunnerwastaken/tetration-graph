"""Prüfung von Fixpunkt und Pentation für alle Basen gegen den Motor (pentx)."""
import mpmath as mp
import lib, motor
mp.mp.dps = 120
for b in ["9/5", "2", "exp(1)", "3", "4", "5", "10"]:
    T = lib.Tet(b)
    m = motor.werte(b, ["pentx(0.5)", "pentx(-0.5)", "pentx(1.5)"], 60)
    own = [T.P(mp.mpf("0.5")), T.P(mp.mpf("-0.5")), T.P(mp.mpf("1.5"))]
    print(f"{b:7s} L={mp.nstr(T.L, 20)} lam={mp.nstr(T.lam, 20)} rho={mp.nstr(T.L+2, 4)}  "
          f"Motor-lib: {[mp.nstr(mp.re(x) - y, 3) for x, y in zip(m, own)]}  P(0.5)={mp.nstr(own[0], 20)}")
motor.schliessen()
