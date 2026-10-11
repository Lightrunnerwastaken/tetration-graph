"""Strang P, Prüfung 3: exakte Identitäten der Hauptkette am echten Problem (f = f_K aus 4000 Koeffizienten, g = exp).

T1  H1 Schritt (1): Lösung (λ,Q) von T_N[Φ((1−λ)f + Q) − eˣ] = 0  ⇒  Q = λ Π^P_N f   (λ = 0.4, 1.0)
T2  H2: P0 = pad(p_N), c_M (T_M 𝒥_{P0})⁻¹ x^M = P0 − Π^{P0}_M f
T3  B §9.5(ii): x1(N+1) − pad x1(N) = −ρ^f_N t_{N+1}
T4  B-Id (beide Gleichungen) mit Zufallspolynomen a, b
T5  Korollar B-Ξ2 / d^(2)-Formel
T6  D8: J_N(f) = M_N (I + K_N)
T7  Lemma 7.1 / C 12.1: y_ζ = ζ^{−N−1} Σ z_m (ζ−a0)^{−m−1}, Abfall E_N ~ ζ^{−N−2}
"""
import random
import sys

import mpmath as mp

from reihen import (add, compose, deriv, div, exp_reihe, horner, inv, kneser_koeff, mul, scal, shift)

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 70
N = 22
K = N + 2
A = kneser_koeff(4000)
a0 = A[0]
EXP = exp_reihe(K)
random.seed(7)


def nrm(v):
    return max(abs(x) for x in v) if len(v) else mp.mpf(0)


def ser(v, K_):
    return (list(v) + [mp.mpf(0)] * (K_ + 1))[:K_ + 1]


# f an 0 und an beliebigen Stellen (Cache)
F0 = ser(A, K + 2)
_shift_cache = {}


def f_at(c, KK):
    key = (mp.nstr(c, 50), KK)
    if key not in _shift_cache:
        _shift_cache[key] = shift(A, c, KK)
    return _shift_cache[key]


def P_data(lam, Q, KK):
    """P = (1−λ)f + Q: Reihe an 0 und Taylor an P(0) (Ordnung KK)."""
    Pz = add(scal(1 - lam, ser(A, KK)), ser(Q, KK), KK)
    c = Pz[0]
    Pc = add(scal(1 - lam, f_at(c, KK)), shift(ser(Q, len(Q) - 1), c, KK), KK)
    return Pz, Pc


def Phi_T(Pz, Pc, KK):
    return compose(Pc, Pz, KK)


def J_mat(Pz, Pc, n):
    """Tangente J_n(P): Spalte j = T_n[P^j + x^j (P'∘P)]."""
    dP = compose(deriv(Pc) + [0], Pz, n)
    cols = []
    pw = ser([1], n)
    for j in range(n + 1):
        col = add(pw, [0] * j + dP[:n + 1 - j], n)
        cols.append(col)
        pw = mul(pw, Pz, n)
    return mp.matrix([[cols[j][k] for j in range(n + 1)] for k in range(n + 1)])


def M_sek(Pz, n):
    """M_P = (f∘P − f∘f)/(P − f) als Reihe bis x^n (P(0) ≠ f(0))."""
    fP = compose(f_at(Pz[0], n + 1), Pz, n)
    ff = compose(f_at(a0, n + 1), ser(A, n), n)
    return div(add(fP, ff, n, -1), add(Pz, ser(A, n), n, -1), n)


def Jsek_mat(Pz, n):
    Mp = M_sek(Pz, n)
    cols = []
    pw = ser([1], n)
    for j in range(n + 1):
        cols.append(add(pw, [0] * j + Mp[:n + 1 - j], n))
        pw = mul(pw, Pz, n)
    return mp.matrix([[cols[j][k] for j in range(n + 1)] for k in range(n + 1)]), Mp


def Pi_sek_f(Pz, n):
    """Π^P_n f = (T_n 𝒥_P|)⁻¹ T_n[f∘P + M_P f]."""
    Jm, Mp = Jsek_mat(Pz, n)
    rhs = add(compose(f_at(Pz[0], n + 1), Pz, n), mul(Mp, ser(A, n), n), n)
    return mp.lu_solve(Jm, mp.matrix(rhs)), Jm


def newton_Q(lam, Q, n, it=30):
    for _ in range(it):
        Pz, Pc = P_data(lam, Q, n + 1)
        F = add(Phi_T(Pz, Pc, n), EXP, n, -1)
        d = mp.lu_solve(J_mat(Pz, Pc, n), mp.matrix(F))
        Q = [Q[k] - d[k] for k in range(n + 1)]
        if nrm(d) < mp.mpf(10) ** (-mp.mp.dps + 8):
            break
    return Q, nrm(F)


print(f"N = {N}, dps = {mp.mp.dps}")
# ---------------- T1
for lam in [mp.mpf("0.4"), mp.mpf(1)]:
    Q, res = newton_Q(lam, scal(lam, ser(A, N)), N)
    Pz, Pc = P_data(lam, Q, N + 1)
    PiF, _ = Pi_sek_f(Pz, N)
    diff = nrm([Q[k] - lam * PiF[k] for k in range(N + 1)])
    print(f"T1 λ = {lam}: Residuum {mp.nstr(res, 3)}, max|Q − λΠ^P f| = {mp.nstr(diff, 3)}, |Q| ~ {mp.nstr(nrm(Q), 3)}")
    if lam == 1:
        pN = Q

# ---------------- T2
Mx = N + 1
P0 = pN + [mp.mpf(0)]
Pz, Pc = P_data(mp.mpf(1), P0, Mx + 1)
cM = Phi_T(Pz, Pc, Mx)[Mx] - 1 / mp.factorial(Mx)
PiF, Jm = Pi_sek_f(Pz, Mx)
eM = mp.matrix([0] * Mx + [1])
lhs = mp.lu_solve(Jm, eM) * cM
diff = nrm([lhs[k] - (P0[k] - PiF[k]) for k in range(Mx + 1)])
print(f"T2 H2-Identität: c_M = {mp.nstr(cM, 4)}, max|c_M 𝒥⁻¹x^M − (P0 − Π^{{P0}} f)| = {mp.nstr(diff, 3)}, "
      f"Grösse {mp.nstr(nrm(lhs), 3)}")
# Pfad-Integral: Homotopie s: 0 → 1 per kleinen Schritten, Vergleich mit Newton-Lösung p_{N+1}
pM, _ = newton_Q(mp.mpf(1), P0, Mx)
print(f"     |p_(N+1) − pad p_N| = {mp.nstr(nrm([pM[k] - P0[k] for k in range(Mx + 1)]), 3)}")


# ---------------- T3
def x1(n):
    """Π_n f = J_n(f)⁻¹ T_n J_∞ f,  J_∞f = eˣ + (f'∘f) f."""
    Fz = ser(A, n + 2)
    Fc = f_at(a0, n + 3)
    dF = compose(deriv(Fc) + [0], Fz, n + 1)
    rhs = add(EXP[:n + 2] + [0] * 3, mul(dF, Fz, n + 1), n + 1)
    Jn = J_mat(Fz, Fc, n)
    x = mp.lu_solve(Jn, mp.matrix(rhs[:n + 1]))
    return [x[k] for k in range(n + 1)], Jn, dF


x1N, JN, dF = x1(N)
x1M, JM, _ = x1(N + 1)
# ρ^f_N = [x^{N+1}] J_∞(x1(N) − f) = [x^{N+1}] (x1∘f + (f'∘f) x1 − eˣ − (f'∘f) f)
Fz = ser(A, N + 2)
x1_at_f = compose(shift(x1N, a0, N + 2), Fz, N + 2)
r = add(add(x1_at_f, mul(dF, x1N, N + 1), N + 1), add(EXP[:N + 2], mul(dF, Fz, N + 1), N + 1), N + 1, -1)
print(f"T3 Kontrolle T_N J_∞(x1 − f) = 0: {mp.nstr(nrm(r[:N + 1]), 3)}")
rho_f = r[N + 1]
tM = mp.lu_solve(JM, mp.matrix([0] * (N + 1) + [1]))
diff = nrm([x1M[k] - (x1N[k] if k <= N else 0) + rho_f * tM[k] for k in range(N + 2)])
print(f"T3 B9.5(ii): ρ^f_N = {mp.nstr(rho_f, 4)}, max|x1(N+1) − pad x1(N) + ρ t_(N+1)| = {mp.nstr(diff, 3)}, "
      f"|x1(N+1) − pad x1(N)| = {mp.nstr(nrm([x1M[k] - (x1N[k] if k <= N else 0) for k in range(N + 2)]), 3)}")

# ---------------- T4  B-Id
KK = 18
fz = ser(A, KK + 2)
fc = f_at(a0, KK + 3)
fp = deriv(fz)
g = compose(deriv(fc) + [0], fz, KK + 1)          # g = f'∘f
fpp_f = compose(deriv(deriv(fc)) + [0, 0], fz, KK)   # f''∘f


def rnd_poly(d):
    return [mp.mpf(random.uniform(-1, 1)) / (k + 1) for k in range(d + 1)]


aP, bP = rnd_poly(7), rnd_poly(7)


def comp_poly_f(p, KK_):
    return compose(shift(p, a0, KK_ + 1), fz, KK_)


def J_inf(p, KK_):
    return add(comp_poly_f(p, KK_), mul(g, p, KK_), KK_)


ra, rb = J_inf(aP, KK + 1), J_inf(bP, KK + 1)
ap_f = comp_poly_f(deriv(aP), KK)
bp_f = comp_poly_f(deriv(bP), KK)
D2 = add(add(mul(mul(fpp_f, aP, KK), bP, KK), mul(ap_f, bP, KK), KK), mul(bp_f, aP, KK), KK)
lhs = mul(fp, D2, KK)
ab = mul(aP, bP, KK + 1)
rhs = add(add(mul(deriv(ra), bP, KK), mul(deriv(rb), aP, KK), KK), deriv(mul(g, ab, KK + 1)), KK, -1)
print(f"T4 B-Id (1): max|f'D²Φ − (r_a'b + r_b'a − (g ab)')| = {mp.nstr(nrm(add(lhs, rhs, KK - 1, -1)), 3)} "
      f"(Grösse {mp.nstr(nrm(lhs), 3)})")
# (2): (g c)'/f' = J_∞[c'/(2f')] + [(g(a r_b + b r_a))' − (r_a r_b)']/(2 g f'),  c = ab
c = ab
lhs2 = div(deriv(mul(g, c, KK + 1)), fp, KK - 1)
v = div(deriv(c), scal(2, fp), KK + 1)                 # v = c'/(2f') an 0
# v an a0: c'(a0+y)/(2 f'(a0+y))
cpa = shift(deriv(c), a0, KK + 1)
fpa = deriv(fc)
va = div(cpa, scal(2, fpa[:KK + 2]), KK + 1)
Jv = add(compose(va, fz, KK), mul(g, v, KK), KK)
term = add(deriv(mul(g, add(mul(aP, rb, KK + 1), mul(bP, ra, KK + 1), KK + 1), KK + 1)), deriv(mul(ra, rb, KK + 1)), KK, -1)
rhs2 = add(Jv, div(term, scal(2, mul(g, fp, KK)), KK), KK)
print(f"T4 B-Id (2): max|…| = {mp.nstr(nrm(add(lhs2, rhs2, KK - 2, -1)), 3)} (Grösse {mp.nstr(nrm(lhs2), 3)})")
# Folge 9.0.1: D²Φ(f)[a,b] = −J_∞[(ab)'/(2f')] + (1/(2f'))[b r_a' − b' r_a + a r_b' − a' r_b − (g'/g)(a r_b + b r_a)] + (r_a r_b)'/(2gf')
gp = deriv(g)
br = add(add(mul(bP, deriv(ra), KK), mul(deriv(bP), ra, KK), KK, -1),
         add(mul(aP, deriv(rb), KK), mul(deriv(aP), rb, KK), KK, -1), KK)
br = add(br, mul(div(gp, g, KK), add(mul(aP, rb, KK), mul(bP, ra, KK), KK), KK), KK, -1)
rhs3 = add(add(scal(-1, Jv), div(br, scal(2, fp), KK), KK), div(deriv(mul(ra, rb, KK + 1)), scal(2, mul(g, fp, KK)), KK), KK)
print(f"T4 Formel 9.0.1: max|D²Φ − rechte Seite| = {mp.nstr(nrm(add(D2, rhs3, KK - 2, -1)), 3)} (Grösse {mp.nstr(nrm(D2), 3)})")

# ---------------- T5  d^(2)
# u1 = x1 − f (an 0 und a0), d2_direct = −J⁻¹ T_N ½ D²Φ(f)[u1,u1],  D²Φ[u,u] = (f''∘f)u² + 2(u'∘f)u
KN = N
u0 = add(ser(x1N, KN + 2), ser(A, KN + 2), KN + 2, -1)
ua = add(shift(x1N, a0, KN + 3), f_at(a0, KN + 3), KN + 3, -1)
fz2 = ser(A, KN + 2)
fc2 = f_at(a0, KN + 4)
g2 = compose(deriv(fc2) + [0], fz2, KN + 2)
fppf2 = compose(deriv(deriv(fc2)) + [0, 0], fz2, KN)
upf = compose(deriv(ua) + [0], fz2, KN)
D2u = add(mul(fppf2, mul(u0, u0, KN), KN), scal(2, mul(upf, u0, KN)), KN)
d2_dir = mp.lu_solve(JN, mp.matrix(scal(-mp.mpf(1) / 2, D2u[:KN + 1])))
# Formel: Π_N[u1 u1'/(2f')] − (N+1) ρ^f u1(0)/(2 f'(0)) t_N
w0 = div(mul(u0, deriv(u0), KN + 1), scal(2, deriv(fz2)), KN)
wa = div(mul(ua, deriv(ua), KN + 1), scal(2, deriv(fc2)[:KN + 3]), KN + 1)
Jw = add(compose(wa, fz2, KN), mul(g2, w0, KN), KN)
Piw = mp.lu_solve(JN, mp.matrix(Jw[:KN + 1]))
tN = mp.lu_solve(JN, mp.matrix([0] * KN + [1]))
ru = r  # J_∞ u1 = r (oben), ρ = r[N+1]
d2_form = [Piw[k] - (KN + 1) * rho_f * u0[0] / (2 * deriv(fz2)[0]) * tN[k] for k in range(KN + 1)]
print(f"T5 d^(2): |direkt| = {mp.nstr(nrm(d2_dir), 3)}, max|direkt − Formel| = "
      f"{mp.nstr(nrm([d2_dir[k] - d2_form[k] for k in range(KN + 1)]), 3)}, Kopfanteil "
      f"{mp.nstr(nrm([(KN + 1) * rho_f * u0[0] / (2 * deriv(fz2)[0]) * tN[k] for k in range(KN + 1)]), 3)}")

# ---------------- T6  D8
n = N
m_ser = dF[:n + 1]                      # m = f'∘f
Mn = mp.matrix([[m_ser[k - j] if k >= j else 0 for j in range(n + 1)] for k in range(n + 1)])
im = inv(m_ser, n)
Minv = mp.matrix([[im[k - j] if k >= j else 0 for j in range(n + 1)] for k in range(n + 1)])
pw = ser([1], n)
Ccols = []
for j in range(n + 1):
    Ccols.append(pw)
    pw = mul(pw, ser(A, n), n)
Cf = mp.matrix([[Ccols[j][k] for j in range(n + 1)] for k in range(n + 1)])
Kn = Minv * Cf
diffD8 = mp.mnorm(JN - Mn * (mp.eye(n + 1) + Kn), 1)
print(f"T6 D8: |J_N − M_N(I+K_N)|_1 = {mp.nstr(diffD8, 3)}, |M_N·Minv − I| = {mp.nstr(mp.mnorm(Mn * Minv - mp.eye(n + 1), 1), 3)}")

# ---------------- T7  Lemma 7.1 / C 12.1
z = mp.mpc("0.3", "0.9")


def Pi_k(zeta):
    # Π_N k_ζ = J_N(f)⁻¹ T_N[1/(ζ − f(x)) + m(x)/(ζ − x)]
    kz0 = [1 / zeta ** (k + 1) for k in range(N + 1)]
    inv_f = inv(add([zeta], ser(A, N), N, -1), N)
    rhs = add(inv_f, mul(dF[:N + 1], kz0, N), N)
    return mp.lu_solve(JN, mp.matrix(rhs))


for zeta in [mp.mpc(3, 1), mp.mpc(12, 5), mp.mpc(50, 20)]:
    Pk = Pi_k(zeta)
    E = horner([Pk[k] for k in range(N + 1)], z) - 1 / (zeta - z)
    print(f"T7 ζ = {mp.nstr(zeta, 4)}: |E_N(ζ,z)|·|ζ|^(N+2) = {mp.nstr(abs(E) * abs(zeta) ** (N + 2), 6)}")
# C 12.1 mit A_f: y_ζ = A_f⁻¹ B_N[k]_{>N} = Π_N k − [k]_{<=N};  Formel ζ^{−N−1} Σ_m z_m (ζ−a0)^{−m−1}, z_m = A⁻¹[q^{N+1} g^m]
q = ser(A, N)
gq = add(q, [a0], N, -1)
qN1 = ser([1], N)
for _ in range(N + 1):
    qN1 = mul(qN1, q, N)
zm = []
gm = ser([1], N)
for m in range(N + 1):
    zm.append(mp.lu_solve(JN, mp.matrix(mul(qN1, gm, N))))
    gm = mul(gm, gq, N)
for zeta in [mp.mpc(1.5, 1.2), mp.mpc(-1.6, 0.3)]:
    Pk = Pi_k(zeta)
    y1 = [Pk[k] - 1 / zeta ** (k + 1) for k in range(N + 1)]
    y2 = [sum(zm[m][k] * (zeta - a0) ** (-m - 1) for m in range(N + 1)) / zeta ** (N + 1) for k in range(N + 1)]
    print(f"T7 C12.1 ζ = {mp.nstr(zeta, 3)}: max|y(direkt) − y(Formel)| = {mp.nstr(nrm([y1[k] - y2[k] for k in range(N + 1)]), 3)}, "
          f"|y| = {mp.nstr(nrm(y1), 3)}")
