"""Gegenproben fuer lib21.sexp: Taylor gegen theta-Darstellung, Funktionalgleichung."""
import numpy as np
import lib21 as S

print("L =", S.L, " P = 2 pi i/L =", S.P_PER)
print("sexp(0), sexp(1)-e, sexp(0.5):", S.sexp(0.0), S.sexp(1.0) - np.e, S.sexp(0.5))
# Taylor gegen theta auf |z|<=1.5, Im z >= 0.25
zs = np.array([0.3 + 0.5j, -0.7 + 0.4j, 1.2 + 0.6j, 0.0 + 1.4j, -1.3 + 0.3j, 0.9 + 0.25j])
tay = np.polyval(S.TAYLOR[::-1], zs)
th = S.F(S.tau(zs))
for z, a, b in zip(zs, tay, th):
    print(f"z={z}: Taylor {a:.15f}  theta {b:.15f}  |diff| {abs(a-b):.2e}")
# Funktionalgleichung
for z in [0.4 + 0.1j, 1.7 + 0.05j, -1.8 + 0.1j, 2.5 + 0.3j, -2.5 + 0.05j]:
    print(f"z={z}: |sexp(z+1)-exp(sexp z)| = {abs(S.sexp(z+1)-np.exp(S.sexp(z))):.2e}")
# Ableitung bei 0 gegen Taylor
print("sexp'(0) num", S.dsexp(0.0), " Taylor", S.TAYLOR[1])
