"""Gewicht der Ebene m in Phi5:  Phi5 - Phi5^(m) = S_m * (c_m - mu4) + O((c_m-mu4)^2),
Vorhersage S_m = (1 + Phi5'(th)) / P_d'(th+m)  (Phi5' grob, h = 1e-4; P_d' exakt numerisch)."""
import mpmath as mp
import lib
mp.mp.dps = 120
Tb, Td = lib.Tet("exp(1)"), lib.Tet("2")
M = 32
mu4 = mp.fsum(lib.phi4(Tb, Td, mp.mpf(j) / M) for j in range(M)) / M
m = 2
for j in range(4, 16):
    th = mp.mpf(j) / 16
    y = Td.P(th)
    U = [y]
    for _ in range(m):
        U.append(Td.T(U[-1]))
    cm = lib.H(Tb, Td, U[m]) - U[m]
    w = U[m] + mu4
    for _ in range(m):
        w = Tb.slog(w)
    phi_m = Tb.spog(w) - th
    phi = lib.phi5(Tb, Td, th)
    h = mp.mpf(10) ** -4
    dphi = (lib.phi5(Tb, Td, th + h) - lib.phi5(Tb, Td, th - h)) / (2 * h)
    e = mp.mpf(10) ** -40
    dP = (Td.P(th + m + e) - Td.P(th + m - e)) / (2 * e)
    S = (1 + dphi) / dP
    print(f"th={float(th):.4f} U_2={mp.nstr(U[m], 6):>9}  (Phi5-Phi5^m)/(c_m-mu4) = {mp.nstr((phi-phi_m)/(cm-mu4), 8):>12}   (1+Phi5')/P_d'(th+m) = {mp.nstr(S, 8):>12}")
