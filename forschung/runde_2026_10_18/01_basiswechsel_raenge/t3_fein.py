"""Feinabtastung von Phi5(e,2) nahe th = 0.3..0.5: Glattheit, Kettenlänge, X-Abhängigkeit."""
import mpmath as mp
import lib

mp.mp.dps = 120
Te, T2 = lib.Tet("exp(1)"), lib.Tet("2")
A = T2.slog(mp.mpf(10) ** 60)
print("A =", mp.nstr(A, 10))
prev = []
for j in range(0, 41):
    th = mp.mpf("0.25") + mp.mpf(j) / 160
    y = T2.P(th)
    U = [y]
    while U[-1] < A:
        U.append(T2.T(U[-1]))
    a = lib.phi5(Te, T2, th, 60)
    b = lib.phi5(Te, T2, th, 40)
    prev.append(a)
    print(f"{float(th):.5f} len={len(U)} Ulast={mp.nstr(U[-1], 5):>10}  Phi5={mp.nstr(a, 20)}  X40-X60={mp.nstr(b-a, 3)}")
d2 = [prev[i + 1] - 2 * prev[i] + prev[i - 1] for i in range(1, len(prev) - 1)]
print("zweite Differenzen:", [mp.nstr(x, 4) for x in d2])
