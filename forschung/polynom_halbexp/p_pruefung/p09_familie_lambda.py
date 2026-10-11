"""Strang P, Prüfung 9: numerische Stütze von (H) an NICHT-polynomialen Punkten des λ-Pfads (bisher nur P = f oder Polynome).
P(λ) = (1−λ)f + Q(λ) mit T_N[Φ(P) − eˣ] = 0 (exakte Lösung, Newton). Gemessen: dist·max_{|z|=ρ}|Π^P_N k_ζ(z)| (wie b25) für
λ = 0, 0.5, 1 und (H-iii): sup|J_N(P)⁻¹xᴺ| / sup|(T_N𝒥_P)⁻¹xᴺ|."""
import sys
import time

import mpmath as mp

from reihen import add, compose, deriv, div, horner, inv, kneser_koeff, mul, scal, shift, fixpunkt_L, exp_reihe

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
A = kneser_koeff(4000)
L = fixpunkt_L()
rho, a0 = abs(L), A[0]
R0 = abs(L - a0)
psi = mp.radians(mp.mpf("87.151205"))
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


Ff = compose(shift(A, a0, N + 1), ser(A, N), N)   # f∘f


CA = shift(A, a0, N + 16)


def M_stabil(Pz, n, J=14):
    """M_P = Σ_{j≥1} (f^{(j)}∘f/j!)(P − f)^{j−1} — ohne Division (die Division ist instabil: P − f hat eine Nullstelle nahe 0)."""
    fz = ser(A, n)
    D = add(Pz, fz, n, -1)
    w = [a0] + add(fz, [a0], n, -1)[1:]
    Mp, Dpow = [mp.mpf(0)] * (n + 1), ser([1], n)
    for j in range(1, J + 1):
        coeffs = [mp.binomial(i + j, j) * CA[i + j] for i in range(n + 1)]
        Mp = add(Mp, mul(compose(coeffs, w, n), Dpow, n), n)
        Dpow = mul(Dpow, D, n)
    return Mp


def Jsek(Pz, n):
    Mp = M_stabil(Pz, n)
    cols, pw = [], ser([1], n)
    for j in range(n + 1):
        cols.append(add(pw, [0] * j + Mp[:n + 1 - j], n))
        pw = mul(pw, Pz, n)
    return mp.matrix([[cols[j][k] for j in range(n + 1)] for k in range(n + 1)]), Mp


def dist(z):
    cands = [abs(z - L), abs(z - mp.conj(L))]
    for (c, r, oc, orr) in [(mp.mpf(0), rho, a0, R0), (a0, R0, mp.mpf(0), rho)]:
        p = c + r * (z - c) / abs(z - c)
        if abs(p - oc) >= orr:
            cands.append(abs(z - p))
    return min(cands)


zs = [rho * mp.expj(2 * mp.pi * k / 360) for k in range(360)]
t0 = time.time()
Q = scal(0, ser(A, N))
for lam in [mp.mpf(0), mp.mpf("0.5"), mp.mpf(1)]:
    Q = scal(lam, ser(A, N)) if lam > 0 else ser([0], N)
    for _ in range(25):
        Pz, Pc = P_data(lam, Q, N + 1)
        F = add(compose(Pc, Pz, N), EXP, N, -1)
        if lam == 0:
            break
        dQ = mp.lu_solve(J_mat(Pz, Pc, N), mp.matrix(F))
        Q = [Q[k] - dQ[k] for k in range(N + 1)]
        if max(abs(dQ[k]) for k in range(N + 1)) < mp.mpf(10) ** -70:
            break
    Pz, Pc = P_data(lam, Q, N + 1)
    res = max(abs(x) for x in add(compose(Pc, Pz, N), EXP, N, -1))
    JS, Mp = Jsek(Pz, N)
    JT = J_mat(Pz, Pc, N)
    vals = []
    for name, basis, richt in [("Ecke", L, mp.expj(psi)), ("links", -rho, mp.mpf(-1))]:
        for t in [mp.mpf("0.5"), mp.mpf(2), mp.mpf(8)]:
            zeta = basis + (t / N) * richt
            kP = inv(add([zeta], Pz[:N + 1], N, -1), N)
            rhs = add(kP, mul(Mp, [1 / zeta ** (k + 1) for k in range(N + 1)], N), N)
            x = mp.lu_solve(JS, mp.matrix(rhs))
            c = [x[k] for k in range(N + 1)]
            vals.append(dist(zeta) * max(abs(horner(c, z)) for z in zs))
    eN = mp.matrix([0] * N + [1])
    tT = mp.lu_solve(JT, eN)
    tS = mp.lu_solve(JS, eN)
    sT = max(abs(horner([tT[k] for k in range(N + 1)], z)) for z in zs)
    sS = max(abs(horner([tS[k] for k in range(N + 1)], z)) for z in zs)
    supPf = max(abs(horner(Pz[:N + 1], z) - horner(A[:600], z)) for z in zs[::6]) if lam > 0 else 0
    print(f"N = {N}, λ = {lam}: Residuum {mp.nstr(res, 3)}, sup_{{|z|=ρ}}|P − f| ≈ {mp.nstr(supPf, 3)} (Reihe 600 Terme)")
    print(f"   dist·max|Π^P k_ζ|  Ecke t=0.5/2/8: {' / '.join(mp.nstr(v, 4) for v in vals[:3])}   links: {' / '.join(mp.nstr(v, 4) for v in vals[3:])}")
    print(f"   (H-iii): sup|J(P)⁻¹xᴺ| = {mp.nstr(sT, 4)}, sup|(T𝒥_P)⁻¹xᴺ| = {mp.nstr(sS, 4)}, Quotient {mp.nstr(sT / sS, 6)}   ({time.time() - t0:.0f} s)")
