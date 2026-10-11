"""Strang P, Prüfung 19: Ms adjungierte Zwei-Punkt-Form (M §4.2) unabhängig an f = f_K (Kneser-Reihe) nachrechnen:
(J_N^T v)_{N−l} = [t^l]( m·W − PP_{a₀}[Φ_W] ),  W(t) = Σ_i v_{N−i} t^i,  Φ_W(y) = W(h(y))·h'(y)·(y/h(y))^{N+1},  m = f'∘f."""
import random
import sys

import mpmath as mp

from reihen import add, compose, deriv, inv, kneser_koeff, mul, shift

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
random.seed(11)
A = kneser_koeff(4000)
a0 = A[0]
N = 14
K = N + 2
# J_N
fz = A[:K + 1]
fc = shift(A, a0, K + 2)
m = compose(deriv(fc) + [0], fz, K)                      # m = f'∘f an 0
cols = []
pw = [mp.mpf(1)] + [mp.mpf(0)] * N
for j in range(N + 1):
    cols.append(add(pw, [0] * j + m[:N + 1 - j], N))
    pw = mul(pw, fz, N)
J = [[cols[j][k] for j in range(N + 1)] for k in range(N + 1)]
v = [mp.mpf(random.uniform(-1, 1)) for _ in range(N + 1)]
JTv = [sum(J[k][j] * v[k] for k in range(N + 1)) for j in range(N + 1)]
W = [v[N - i] for i in range(N + 1)]                    # W(t) = Σ v_{N−i} t^i
# h = f⁻¹ bei a0: Reihe in y − a0 (Umkehrung von f(x) − a0 = Σ_{k≥1} a_k x^k), Ordnung K
Mk = K + 3
fa = [mp.mpf(0)] + A[1:Mk + 1]
# Lagrange/Newton-Umkehrung: h(w) mit f(h(w)) − a0 = w
h = [mp.mpf(0), 1 / fa[1]] + [mp.mpf(0)] * (Mk - 1)
for _ in range(40):
    comp = compose(fa, h, Mk)                            # f(h) − a0 als Reihe in w
    err = add(comp, [0, 1], Mk, -1)
    dfh = compose(deriv(fa) + [0], h, Mk)
    h = add(h, mul(err, inv(dfh, Mk), Mk), Mk, -1)
# Φ_W(y) mit y = a0 + w: W(h(w))·h'(w)·((a0+w)/h(w))^{N+1}; h(w) = w·g(w), g(0) = 1/f'(0)
g = h[1:] + [mp.mpf(0)]
Wh = compose(W + [0] * 3, h, Mk)
hp = deriv(h)
y = [a0, mp.mpf(1)]
q = mul(y, inv(g, Mk), Mk)                               # (a0 + w)/g(w)
qN = [mp.mpf(1)] + [mp.mpf(0)] * Mk
for _ in range(N + 1):
    qN = mul(qN, q, Mk)
ser = mul(mul(Wh, hp, Mk), qN, Mk)                       # = Φ_W · w^{N+1}  (Potenzreihe in w)
# Hauptteil: Σ_{n=1}^{N+1} ser[N+1−n] w^{−n}; Entwicklung in t (w = t − a0): w^{−n} = (−a0)^{−n} Σ_l binom(n+l−1,l) (t/a0)^l
PP = [mp.mpf(0)] * (N + 1)
for n in range(1, N + 2):
    c = ser[N + 1 - n]
    for l in range(N + 1):
        PP[l] += c * (-a0) ** (-n) * mp.binomial(n + l - 1, l) * a0 ** (-l)
mW = mul(m, W, N)
rhs = [mW[l] - PP[l] for l in range(N + 1)]
lhs = [JTv[N - l] for l in range(N + 1)]
print(f"N = {N}: max |(J^T v)_{{N−l}} − [t^l](mW − PP)| = {mp.nstr(max(abs(lhs[l] - rhs[l]) for l in range(N + 1)), 3)}, "
      f"max|lhs| = {mp.nstr(max(abs(x) for x in lhs), 3)}")
