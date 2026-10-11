"""Strang P, Prüfung 9b: Kontrolle von p09 — M_P = (f∘P − f∘f)/(P − f) per Division gegen die stabile Reihe
M_P = Σ_{j≥1} (f^{(j)}∘f / j!) (P − f)^{j−1}; danach Kopfantwort sup|(T_N𝒥_P)⁻¹xᴺ| mit beiden Formen."""
import sys

import mpmath as mp

from reihen import add, compose, deriv, div, horner, kneser_koeff, mul, scal, shift, fixpunkt_L, exp_reihe

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
A = kneser_koeff(4000)
L = fixpunkt_L()
rho, a0 = abs(L), A[0]
EXP = exp_reihe(N + 2)


def ser(v, K):
    return (list(v) + [mp.mpf(0)] * (K + 1))[:K + 1]


def P_data(lam, Q, K):
    Pz = add(scal(1 - lam, ser(A, K)), ser(Q, K), K)
    c = Pz[0]
    Pc = add(scal(1 - lam, shift(A, c, K)), shift(ser(Q, len(Q) - 1), c, K), K)
    return Pz, Pc


def J_mat(Pz, Pc, n):
    dP = compose(deriv(Pc) + [0], Pz, n)
    cols, pw = [], ser([1], n)
    for j in range(n + 1):
        cols.append(add(pw, [0] * j + dP[:n + 1 - j], n))
        pw = mul(pw, Pz, n)
    return mp.matrix([[cols[j][k] for j in range(n + 1)] for k in range(n + 1)])


lam = mp.mpf("0.5")
Q = scal(lam, ser(A, N))
for _ in range(25):
    Pz, Pc = P_data(lam, Q, N + 1)
    F = add(compose(Pc, Pz, N), EXP, N, -1)
    dQ = mp.lu_solve(J_mat(Pz, Pc, N), mp.matrix(F))
    Q = [Q[k] - dQ[k] for k in range(N + 1)]
    if max(abs(dQ[k]) for k in range(N + 1)) < mp.mpf(10) ** -70:
        break
Pz, Pc = P_data(lam, Q, N + 1)
fz = ser(A, N)
D = add(Pz, fz, N, -1)
# Division
fP = compose(shift(A, Pz[0], N + 1), Pz, N)
ff = compose(shift(A, a0, N + 1), fz, N)
M_div = div(add(fP, ff, N, -1), D, N)
# stabile Reihe: Taylor von f an a0 bis Ordnung J+N
J = 12
c = shift(A, a0, N + J + 2)
M_st = [mp.mpf(0)] * (N + 1)
Dpow = ser([1], N)
w = add(fz, [a0], N, -1)                  # f(x) − a0
for j in range(1, J + 1):
    # f^{(j)}(a0 + w)/j! = Σ_i binom(i+j, j) c_{i+j} w^i
    coeffs = [mp.binomial(i + j, j) * c[i + j] for i in range(N + 1)]
    fj_f = compose(coeffs, [a0] + w[1:], N)
    M_st = add(M_st, mul(fj_f, Dpow, N), N)
    Dpow = mul(Dpow, D, N)
print("max|D_k|ρ^k (P − f) =", mp.nstr(max(abs(D[k]) * rho ** k for k in range(N + 1)), 3), "  D(0) =", mp.nstr(D[0], 3))
rel = max(abs(M_div[k] - M_st[k]) * rho ** k for k in range(N + 1))
print("max_k |M_div − M_stabil|·ρ^k =", mp.nstr(rel, 3), "  max|M|ρ^k =", mp.nstr(max(abs(M_st[k]) * rho ** k for k in range(N + 1)), 3))
mff = compose(deriv(shift(A, a0, N + 2)) + [0], fz, N)
print("max_k |M_stabil − f'∘f|·ρ^k =", mp.nstr(max(abs(M_st[k] - mff[k]) * rho ** k for k in range(N + 1)), 3))
zs = [rho * mp.expj(2 * mp.pi * k / 360) for k in range(360)]
for name, Mp in [("Division", M_div), ("stabil", M_st), ("f'∘f", mff)]:
    cols, pw = [], ser([1], N)
    for jj in range(N + 1):
        cols.append(add(pw, [0] * jj + Mp[:N + 1 - jj], N))
        pw = mul(pw, Pz, N)
    JS = mp.matrix([[cols[jj][k] for jj in range(N + 1)] for k in range(N + 1)])
    t = mp.lu_solve(JS, mp.matrix([0] * N + [1]))
    print(f"  {name:9s}: sup|(T_N𝒥)⁻¹xᴺ| = {mp.nstr(max(abs(horner([t[k] for k in range(N + 1)], z)) for z in zs), 5)}")
