"""Konjugator H nahe dem Fixpunkt: G^(x) = spog_b(H(P_d(x))) = x + Phi5({x}) auch für x -> -oo,
H(P_d(x)) -> L_b; Schröder-Form: sigma_b(H(y)) = C sigma_d(y)^kappa exp(ln lam_b * P5(spog_d y)), kappa = ln lam_b / ln lam_d."""
import mpmath as mp
import lib
mp.mp.dps = 120
Tb, Td = lib.Tet("exp(1)"), lib.Tet("2")
p = lib.phi5(Tb, Td, mp.mpf("0.3"))
for x in ["0.3", "-1.7", "-3.7", "-7.7", "-15.7"]:
    x = mp.mpf(x)
    y = Td.P(x)
    Hy = lib.H(Tb, Td, y)
    print(f"x={float(x):6.1f}  P_d(x)-L_d={mp.nstr(y - Td.L, 5):>12}  H-L_b={mp.nstr(Hy - Tb.L, 5):>12}  "
          f"spog_b(H(P_d x)) - x - Phi5(0.3) = {mp.nstr(Tb.spog(Hy) - x - p, 3)}   G(x)-x={mp.nstr(Tb.spog(y) - x if y > Tb.L else mp.nan, 8)}")
print("kappa = ln lam_b/ln lam_d =", mp.nstr(mp.log(Tb.lam) / mp.log(Td.lam), 20), "  spog_b(L_d) =", mp.nstr(Tb.spog(Td.L), 20))
