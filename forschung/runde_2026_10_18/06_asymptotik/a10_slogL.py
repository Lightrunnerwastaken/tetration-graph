"""A10: slog nahe L (Umkehrung von A9).  Für w = sexp(z), Im z gross:
 z = slog(w) erfüllt  z + theta(z) - a0 = (log(sigma(w)/kappa) + 2 pi i j)/L  für ein ganzes j  (Zweigwahl = Blatt von slog).
 Daraus: slog(w) = l/L - a1 (sigma(w)/kappa)^gamma e^{2 pi i j gamma ...}... ; Leitordnung slog(w) = l/L + O(|w-L|^{Re gamma}).
Prüfung: z vorgeben, w = sexp(z) per Motor, z aus der Formel zurückgewinnen.
"""
import mpmath as mp
from motor import sexp
from lib3_kopie import Basis, lade_theta
mp.mp.dps = 50
B = Basis(60); L = B.L; a = lade_theta(); kap = mp.exp(L * a[0])
gam = 2j * mp.pi / L
print("gamma = 2 pi i/L =", mp.nstr(gam, 15))
def inv(l):
    z = l / L
    for _ in range(80):
        q = mp.expjpi(2 * z); z = l / L - mp.fsum(a[n] * q ** n for n in range(1, 40))
    return z
for z in [mp.mpc(0.3, 1), mp.mpc(0.7, 1.5), mp.mpc(-0.4, 2), mp.mpc(0.1, 3)]:
    w = sexp([z], 40)[0]
    l0 = mp.log(B.sigma(w)[0] / kap)
    j = int(mp.nint(mp.re((L * z - l0) / (2j * mp.pi))))
    l = l0 + 2j * mp.pi * j
    z1 = l / L
    z2 = z1 - a[1] * mp.expjpi(2 * z1)
    z3 = inv(l)
    print(f"z={mp.nstr(z,4)} |w-L|={mp.nstr(abs(w-L),3)} j={j}: |z-l/L|={mp.nstr(abs(z-z1),3)}  |z-(l/L-a1 q)|={mp.nstr(abs(z-z2),3)}  |z-voll|={mp.nstr(abs(z-z3),3)}")
