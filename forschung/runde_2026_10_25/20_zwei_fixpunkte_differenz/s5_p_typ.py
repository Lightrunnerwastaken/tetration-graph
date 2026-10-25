"""s5: Singularitätstyp von p bei θ* = α_a(1).
Modell (Herleitung im BERICHT): p(θ) ≈ α_A^re(1/(κ(θ-θ*))) + const, α_A^re = reelle Abel-Funktion von b^x auf (A, ∞)
(reguläre Abel-Funktion am abstossenden Fixpunkt, rechts von A). Daraus
    k |c_k(p)| e^{S_a k} ≈ K · α_A^re'(X0 + ln k / ln b) / ln b
mit zwei freien Konstanten X0, K. Test: X0 aus k=1 und k=kmax bestimmen, alle übrigen k vorhersagen.
Gegenmodelle: K k^γ (Potenz, 2 Parameter) und K/(ln k + C) (2 Parameter).
Aufruf: python s5_p_typ.py t dps N K kmax
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from zf import Paar

t, dps, N, K, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
klo = int(sys.argv[6]) if len(sys.argv) > 6 else 1
mp.mp.dps = dps
P = Paar(mp.mpf(t), K=K)
Sa = P.S_a()
cs = P.fourier(P.p, N)
y = {k: k * abs(cs[k]) * mp.exp(Sa * k) for k in range(1, kmax + 1)}


def alpha_re(X):
    n = 0
    while abs(X - P.A) > P.eps:
        X = P.g(X)
        n += 1
    u = P.sigma(X, P.A, P.qA, P.dqA)
    return n + mp.log(u) / P.LlA


def dalpha(X):
    return mp.diff(alpha_re, X)


mp.mp.dps = 30
lb = P.c


def model(X0, k):
    return dalpha(X0 + mp.log(k) / lb)


# X0 aus Verhältnis k=1 zu k=kmax
target = y[kmax] / y[klo]
f = lambda X0: model(X0, kmax) / model(X0, klo) - target
# Startwert: Gitter
best = None
for j in range(1, 200):
    X0 = P.A + mp.mpf(j) * P.A / 50
    try:
        v = f(X0)
    except Exception:
        continue
    if best is None or abs(v) < abs(best[1]):
        best = (X0, v)
X0 = mp.findroot(f, best[0])
Kc = y[klo] / model(X0, klo)
print(f"t={t} b={mp.nstr(P.b, 12)} A={mp.nstr(P.A, 12)} ln b={mp.nstr(lb, 10)}  X0={mp.nstr(X0, 12)}  K={mp.nstr(Kc * lb, 10)}")
# Gegenmodell Potenz und 1/(ln k + C) mit denselben zwei Stützstellen
gam = mp.log(target) / mp.log(mp.mpf(kmax) / klo)
Cl = mp.findroot(lambda C: (mp.log(klo) + C) / (mp.log(kmax) + C) - target, 1)
print(" k   y_k=k|c_k|e^{S_a k}   slog-Modell/y-1   Potenz k^%s /y-1   1/(ln k+C) /y-1" % mp.nstr(gam, 5))
for k in range(1, kmax + 1):
    m1 = Kc * model(X0, k)
    m2 = y[klo] * (mp.mpf(k) / klo) ** gam
    m3 = y[klo] * (mp.log(klo) + Cl) / (mp.log(k) + Cl)
    print(f"{k:3d} {mp.nstr(y[k], 12):>16} {mp.nstr(m1 / y[k] - 1, 3):>12} {mp.nstr(m2 / y[k] - 1, 3):>12} {mp.nstr(m3 / y[k] - 1, 3):>12}")
