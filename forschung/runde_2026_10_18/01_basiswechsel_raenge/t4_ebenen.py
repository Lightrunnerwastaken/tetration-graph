"""Ebenen-Zerlegung der Rang-5-Phase.

Phi5(th) = spog_b(H(y)) - th, y = P_d(th), U_m = T_d^m(y) = P_d(th+m),
H(y) = slog_b^m(U_m + c_m),  c_m = H(U_m) - U_m  (exakt).
Ersetzt man c_m durch das Mittel mu4, entsteht Phi5^(m). Vorhersage:
Phi5 - Phi5^(m) ~ S_m * P4(U_m),  S_m = spog_b'(H(y)) * (slog_b^m)'(U_m)  (~ 1/P_d'(th+m)).
"""
import sys
import mpmath as mp
import lib

mp.mp.dps = 120
bb, dd = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("exp(1)", "2")
Tb, Td = lib.Tet(bb), lib.Tet(dd)
M = 32
mu4 = mp.fsum(lib.phi4(Tb, Td, mp.mpf(j) / M) for j in range(M)) / M
print("mu4", mp.nstr(mu4, 30))


def zerlegung(th, m):
    y = Td.P(th)
    U = [y]
    for _ in range(m):
        U.append(Td.T(U[-1]))
    Hy = lib.H(Tb, Td, y)
    cm = lib.H(Tb, Td, U[m]) - U[m]
    # slog_b^m an U_m + mu4, mit Ableitung an U_m
    w, ab = U[m] + mu4, mp.mpf(1)
    w2 = U[m]
    for _ in range(m):
        w = Tb.slog(w)
        s, ds = Tb.slog(w2, True)
        ab *= ds
        w2 = s
    phi_m = Tb.spog(w) - th
    _, dsp = Tb.spog(Hy, True)
    phi = Tb.spog(Hy) - th
    P4 = lib.phi4(Tb, Td, U[m] - mp.floor(U[m])) - mu4
    return phi, phi_m, dsp * ab, P4, cm, U[m]


m = int(sys.argv[3]) if len(sys.argv) > 3 else 2
print(f"Ebene m = {m}")
for j in range(0, 17):
    th = mp.mpf(j) / 16
    phi, phi_m, S, P4, cm, Um = zerlegung(th, m)
    print(f"th={float(th):.4f} U_m={mp.nstr(Um, 6):>9} c_m-mu4={mp.nstr(cm-mu4, 6):>12} P4(U_m)={mp.nstr(P4, 6):>12} "
          f"Phi5-Phi5^m={mp.nstr(phi-phi_m, 8):>14}  S*(c_m-mu4)={mp.nstr(S*(cm-mu4), 8):>14}  S={mp.nstr(S, 5)}")
