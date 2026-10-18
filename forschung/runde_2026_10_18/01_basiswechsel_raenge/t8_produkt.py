"""Produktformel H'(y) = prod_j F'(U_j) (1+o(1)), F = slog_b o T_d, U_j = T_d^j(y):
Quotienten H_{n+1}'(y)/H_n'(y) gegen F'(U_n) = 1 + Phi4'(U_n) + o(1); dazu Differenzenquotienten von Phi5 über viele Skalen."""
import mpmath as mp
import lib, lib6
mp.mp.dps = 120
Tb, Td = lib.Tet("exp(1)"), lib.Tet("2")
eps = mp.mpf(10) ** -45
def abl(f, x):
    return (f(x + eps) - f(x - eps)) / (2 * eps)
for th in ["0.05", "0.3", "0.5", "0.7"]:
    th = mp.mpf(th)
    y = Td.P(th)
    U = [y]
    for _ in range(3):
        U.append(Td.T(U[-1]) if U[-1] < 6 else mp.inf)
    Hp = [abl(lambda z: lib6.Hn(Tb, Td, z, n), y) for n in range(1, 5)]
    print(f"th={float(th)}  U = {[mp.nstr(u, 6) for u in U]}")
    for n in range(1, 4):
        if U[n] == mp.inf or U[n] > mp.mpf(10) ** 50:
            break
        Fp = abl(lambda z: lib.R(Tb, Td, z), U[n])
        p4p = abl(lambda z: lib.phi4(Tb, Td, z), U[n] - mp.floor(U[n]))
        print(f"   n={n}: H_{n+1}'/H_{n}' - 1 = {mp.nstr(Hp[n]/Hp[n-1] - 1, 8):>15}   F'(U_n)-1 = {mp.nstr(Fp - 1, 8):>15}   Phi4'(U_n) = {mp.nstr(p4p, 8):>15}")
# Differenzenquotienten von Phi5 bei th0 = 0.05 über die Skalen h = 2^-k
th0 = mp.mpf("0.05")
print("Differenzenquotienten (Phi5(th0+h)-Phi5(th0-h))/(2h), th0 = 0.05:")
prev = None
for k in list(range(2, 40, 2)) + [50, 60, 80]:
    h = mp.mpf(2) ** -k
    D = (lib.phi5(Tb, Td, th0 + h) - lib.phi5(Tb, Td, th0 - h)) / (2 * h)
    print(f"  h=2^-{k:<3d} D = {mp.nstr(D, 15)}" + (f"   Änderung {mp.nstr(D - prev, 4)}" if prev is not None else ""))
    prev = D
