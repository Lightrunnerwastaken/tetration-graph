"""E8b: (i) Zweigpunkt-Kandidat w aus E8 bei dps 60 nachrechnen; (ii) unendliche Familie von Lösungen sexp(w) = 2 pi i
über die Kette v = log^[k](2 pi i) -> L und die Theta-Umkehr z_j = inv(log(sigma(v)/kappa) + 2 pi i j), w_j = z_j + k.
"""
import mpmath as mp
from motor import sexp
from lib3_kopie import Basis, lade_theta
mp.mp.dps = 70
w0 = mp.mpc("2.004312482594688003944313", "0.7094135408694036266790701")
def newton(w, c, dps):
    for _ in range(30):
        h = mp.mpf(10) ** (-dps // 3)
        a, b, f0 = sexp([w + h, w - h, w], dps)
        st = (f0 - c) / ((a - b) / (2 * h)); w -= st
        if abs(st) < mp.mpf(10) ** (-dps + 8): break
    return w
c = 2j * mp.pi
wA = newton(w0, c, 60); wB = newton(w0, c, 45)
print("w (dps60) =", mp.nstr(wA, 40)); print("|w60-w45| =", mp.nstr(abs(wA - wB), 3))
print("sexp(w-1) =", mp.nstr(sexp([wA - 1], 60)[0], 30), " log(2 pi i) =", mp.nstr(mp.log(c), 30))
print("sexp'(w) ~", mp.nstr((sexp([wA + 1e-20], 60)[0] - c) / 1e-20, 15))
mp.mp.dps = 50
B = Basis(60); L = B.L; a = lade_theta(); kap = mp.exp(L * a[0])
def inv(l):
    z = l / L
    for _ in range(100):
        q = mp.expjpi(2 * z); z = l / L - mp.fsum(a[n] * q ** n for n in range(1, 40))
    return z
k = 12
v = c
for _ in range(k): v = mp.log(v)
print("v = log^[12](2 pi i) =", mp.nstr(v, 15), " |v-L| =", mp.nstr(abs(v - L), 3))
l0 = mp.log(B.sigma(v)[0] / kap)
for j in range(-1, 4):
    z = inv(l0 + 2j * mp.pi * j)
    if mp.im(z) < 0.8: print("j", j, "Im z zu klein", mp.nstr(z, 8)); continue
    sv = sexp([z], 40)[0]
    # vorwärts: sexp(z + k) = exp^k(sexp z)
    t = sv
    for _ in range(k): t = mp.exp(t)
    # auf möglichst kleines Re zurückschieben: w = z + k, aber sexp(w - m) = log^[m] Kette (nur Hauptzweig gültig, wenn Kette stetig)
    print(f"j={j}: z={mp.nstr(z,12)}  |sexp(z)-v|={mp.nstr(abs(sv-v),3)}  w=z+{k}={mp.nstr(z+k,12)}  exp^k(sexp z)-2pi i={mp.nstr(abs(t-c),3)}")
