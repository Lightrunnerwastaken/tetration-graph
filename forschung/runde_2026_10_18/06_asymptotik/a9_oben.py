"""A9: Asymptotik oben. Kneser: sexp(z) = Psi(kappa e^{Lz} exp(L(theta(z)-a0))), theta = sum a_n e^{2 pi i n z}.
Volle Doppelreihe: sexp(z) = L + sum_{k>=1} sum_{m>=0} C_{k,m} e^{kLz} e^{2 pi i m z},  C_{k,0} = p_k kappa^k  (Poincare-/Schroeder-Umkehr-Koeffizienten).
Prüfung: Fehler von (i) L+kappa e^{Lz}, (ii) Psi(kappa e^{Lz}) (alle m=0-Terme), (iii) + erste Fourier-Korrektur, gegen Motor (dps 40) und Theta-Reihe (dps 44).
"""
import mpmath as mp
from motor import sexp
from lib3_kopie import Basis, lade_theta
mp.mp.dps = 50
B = Basis(60); L = B.L
a = lade_theta(); kap = mp.exp(L * a[0])
print("p_k (Psi-Koeffizienten), C_{k,0} = p_k kappa^k:")
for k in range(1, 7): print(k, mp.nstr(B.P[k], 15), mp.nstr(B.P[k] * kap ** k, 15))
x = mp.mpf("0.3")
for y in [0.5, 1, 1.5, 2, 3, 4]:
    z = mp.mpc(x, y)
    m = sexp([z], 40)[0]
    U = kap * mp.exp(L * z); q = mp.expjpi(2 * z)
    e1 = L + U
    e2 = B.Psi(U)
    w, d = B.Psi(U, True)
    e3 = w + d * U * L * a[1] * q
    print(f"y={y}: |m-(L+kU)|={mp.nstr(abs(m-e1),3)}  |m-Psi(kU)|={mp.nstr(abs(m-e2),3)}  |m-[..+Fourier1]|={mp.nstr(abs(m-e3),3)}   |U|={mp.nstr(abs(U),3)} |q|={mp.nstr(abs(q),3)}")
