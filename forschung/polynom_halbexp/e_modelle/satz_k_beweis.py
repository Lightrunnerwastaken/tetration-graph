"""Satz K, explizites Modell — rechnergestützter Komplettbeweis (arb-Intervallrechnung).

Modell: f(x) = 1/(3 − x), g = f∘f = (3 − x)/(8 − 3x). Behauptung: Der von N = 1 aus fortgesetzte Zweig p_N (exakt gedachte
Newton-Iteration vom aufgefüllten Vorgänger p_{N−1}, a_N = 0) existiert für ALLE N, und ‖p_N − T_N f‖_1 ≤ r_−(N) = O(6^{−N});
insbesondere p_N → f gleichmässig auf |x| ≤ 1.

Norm ‖δ‖ = Σ|δ_k| (W_1). F_N(p) = T_N(p∘p − g): P_N → P_N. Zulässige Menge U = {p : ‖p − f‖ ≤ ρ_U}, ρ_U = R + 1/18, dort
‖p‖ ≤ r2 = 1/2 + ρ_U < 1 und F_N' ist γ-Lipschitz, ‖F_N'‖ ≤ Lip_F (N-unabhängig). Kantorovich (Ortega–Rheinboldt 12.6.2):
β ≥ ‖F'(x0)⁻¹‖, η ≥ ‖F'(x0)⁻¹F(x0)‖, h = βγη ≤ 1/2 ⇒ Lösung in B̄(x0, r_−), eindeutig in B(x0, r_+) ∩ (Lipschitz-Gebiet),
Newton von x0 konvergiert dorthin; r_∓ = (1 ∓ √(1−2h))/(βγ).
Teil 1: ‖A_N⁻¹‖ (Linearisierung an f) für N ≥ n0 über Blockzerlegung (Satz K (c)).  Teil 2: Kantorovich um q_N = T_N f für
alle N (N ≤ n0+1 direkt, N ≥ n0 über Teil 1; Schranken monoton in N).  Teil 3: der Zweig trifft genau diese Lösungen.
"""
import sys
from flint import arb, arb_mat, arb_poly, ctx
from kern import F_J, zweig

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 400
DREI = arb(3)
n0 = int(sys.argv[1]) if len(sys.argv) > 1 else 16
KMAX = 400
R = arb("0.05")
rhoU = R + arb(1) / 18
r2 = arb("0.5") + rhoU
up = lambda x: x.upper()
amax = lambda xs: arb(max(up(x) for x in xs))
amin = lambda xs: arb(min(up(x) for x in xs))


def fpoly(L):
    return arb_poly([DREI ** (-(k + 1)) for k in range(L)])


def norm1(P, L):
    return sum((abs(P[k]) for k in range(L)), arb(0))


def g_koef(N):      # g = (3 − x)·(1/8)Σ(3x/8)^k
    return [arb(3) / 8] + [arb(3) / 8 * (arb(3) / 8) ** k - (arb(1) / 8) * (arb(3) / 8) ** (k - 1) for k in range(1, N + 1)]


# ---------------- Teil 1 ----------------
f = fpoly(KMAX)
e3f2 = ((DREI - f) * (DREI - f)).truncate(KMAX)          # (3 − f)² = 1/m, m = f'∘f
tail15 = arb("1.5") ** (-KMAX)
w15 = (DREI + arb(2) / 3) ** 2                           # ‖(3 − f)²‖_{1.5}; ‖f^l‖_{1.5} = (2/3)^l
mu1 = norm1(e3f2, KMAX) + tail15 * w15
fj = [arb_poly([1])]
for j in range(1, n0 + 1):
    fj.append((fj[-1] * f).truncate(KMAX))
Kcol = [(fj[j] * e3f2).truncate(KMAX) for j in range(n0 + 1)]      # K x^j = f^j (3 − f)²
IPKP = arb_mat(n0 + 1, n0 + 1)
for i in range(n0 + 1):
    for j in range(n0 + 1):
        IPKP[i, j] = Kcol[j][i] + (1 if i == j else 0)
V = IPKP.inv()
spalten = [arb(1)]
for j in range(n0 + 1):
    s, tl = arb_poly([0]), arb(0)
    for l in range(n0 + 1):
        s += V[l, j] * Kcol[l]
        tl += abs(V[l, j]) * tail15 * w15 * (arb(2) / 3) ** l
    spalten.append(norm1(arb_poly([0] * j + [1]) - s, KMAX) + tl)
B_low = amax(spalten)                                     # ‖(I + XP)⁻¹‖, XP = T_N K P
eps = mu1 * arb(2) ** (-(n0 + 1))                         # ‖X Q‖ ≤ sup_{k>n0} 2^{−k} μ1
assert up(B_low * eps) < 1
Ainv = mu1 * B_low / (1 - B_low * eps)
print(f"Teil 1: μ1 = ‖(3−f)²‖ ≤ {up(mu1).str(7)}, B_low ≤ {up(B_low).str(7)}, ε ≤ {up(eps).str(4)} ⇒ ‖A_N⁻¹‖ ≤ {up(Ainv).str(7)} für alle N ≥ {n0}")

# Lipschitz-Konstanten auf U
S1 = amax(k * r2 ** (k - 1) for k in range(1, 400))
S2b = amax(k * (k - 1) * r2 ** (k - 2) for k in range(2, 600))
gamma = 2 * S1 + arb(2) / (DREI - r2) ** 3 + S2b * rhoU
Lip_F = 1 + arb(1) / (DREI - r2) ** 2 + S1 * rhoU
print(f"auf U (ρ_U = {up(rhoU).str(4)}, r2 = {up(r2).str(4)}): γ ≤ {up(gamma).str(5)}, ‖F'‖ ≤ {up(Lip_F).str(5)}")


def qf(N):
    return DREI ** (-(N + 1)) / 2                         # ‖q_N − f‖


def resid(N):
    return arb(2) / 5 * arb(6) ** (-(N + 1))              # ‖F_N(q_N)‖ ≤ Σ_{k>N} 3^{−k−1}(1/2)^k


def beta_gross(N):
    eN = qf(N) + sum((k * DREI ** (-(k + 1)) * arb(2) ** (-(k - 1)) for k in range(N + 1, N + 400)), arb(0)) \
         + arb(2) ** (-(N + 398)) + arb(2) / arb("2.5") ** 3 * qf(N)
    assert up(Ainv * eN) < 1
    return Ainv / (1 - Ainv * eN)


def kanto(beta, eta):
    h = beta * gamma * eta
    assert up(h) <= arb("0.5"), f"h = {h}"
    w = (1 - 2 * h).sqrt()
    return h, (1 - w) / (beta * gamma), (1 + w) / (beta * gamma)


def direkt(a, N):
    """β, η, ‖F‖ an der Stelle a (Liste arb) für System N."""
    F, J = F_J(a, g_koef(N), N)
    Ji = J.inv()
    b = amax(sum((abs(Ji[i, j]) for i in range(N + 1)), arb(0)) for j in range(N + 1))
    d = Ji * F
    return b, sum((abs(d[i, 0]) for i in range(N + 1)), arb(0)), sum((abs(F[i, 0]) for i in range(N + 1)), arb(0))


# ---------------- Teil 2: Kantorovich um q_N ----------------
K2 = {}
for N in range(1, n0 + 2):
    q = [DREI ** (-(k + 1)) for k in range(N + 1)]
    b, eta, _ = direkt(q, N)
    if N >= n0:
        b = amin([b, beta_gross(N)])
    h, rm, rp = kanto(b, eta)
    assert up(rm) < up(R)
    K2[N] = (b, eta, h, rm, amin([rp, R]))
bG = beta_gross(n0)
hG, rmG, rpG = kanto(bG, bG * resid(n0))
assert up(rmG) < up(R)
print(f"Teil 2: für jedes N gibt es genau eine Lösung p_N* in B(q_N, min(r_+, R)), ‖p_N* − q_N‖ ≤ r_−(N).")
for N in (1, 2, 3, 5, 10, n0):
    b, eta, h, rm, rp = K2[N]
    print(f"   N={N:2d}: β ≤ {up(b).str(5)}, η ≤ {up(eta).str(4)}, h ≤ {up(h).str(4)}, r_− ≤ {up(rm).str(4)}, r_+ ≥ {rp.lower().str(4)}")
print(f"   N ≥ {n0} (monoton): β ≤ {up(bG).str(5)}, h ≤ {up(hG).str(4)}, r_−(N) ≤ {up(bG * arb('0.8')).str(4)}·6^(−N−1), r_+ ≥ {rpG.lower().str(4)}")

# ---------------- Teil 3: der Zweig ----------------
# N = 1: p_1 = 1/3 + x/8 ist exakte Lösung (a_1 = √g_1 = 1/8, a_0 = g_0/(1 + a_1) = 1/3); ‖p_1 − q_1‖ = 1/72.
assert up(abs(arb(1) / 8 - arb(1) / 9)) < K2[1][4].lower()
num = zweig(lambda N: g_koef(N), n0 + 1, bits=400)      # numerische Zentren p̃_N
cent, epsN = {}, {}
for N in range(1, n0 + 1):
    c = [arb(str(x)) for x in num[N]]
    c = [arb(x.mid()) for x in c]
    b, eta, _ = direkt(c, N)
    h, rm, rp = kanto(b, eta)
    dq = sum((abs(c[k] - DREI ** (-(k + 1))) for k in range(N + 1)), arb(0))
    assert up(dq + K2[N][3]) < rp.lower() and up(dq + K2[N][3]) < up(R)   # p_N* liegt in der Eindeutigkeitskugel um p̃_N
    cent[N], epsN[N] = c, rm                                               # ⇒ ‖p_N* − p̃_N‖ ≤ r_−(p̃_N)
schritte = {}
for N in range(1, n0):
    # Box um das aufgefüllte p_N* (enthält es sicher), dann Intervall-Newton-Schritte für System N+1, bis Kantorovich greift
    X = [arb(x.mid(), up(epsN[N])) for x in cent[N]] + [arb(0)]
    gk1 = g_koef(N + 1)
    for it in range(8):
        c = [arb(x.mid()) for x in X]
        e = sum((arb(x.rad()) for x in X), arb(0))                 # Box ⊂ l¹-Kugel B(c, e)
        b_c, eta_c, _ = direkt(c, N + 1)
        t = b_c * gamma * e
        ok = up(t) < 1
        if ok:
            b0 = b_c / (1 - t)
            eta0 = eta_c / (1 - t) + b0 * Lip_F * e
            h0 = b0 * gamma * eta0
            ok = up(h0) <= arb("0.5")
        if ok:
            h0, rm0, _ = kanto(b0, eta0)
            dz = sum((abs(c[k] - DREI ** (-(k + 1))) for k in range(N + 2)), arb(0)) + e + rm0
            ok = up(dz) < K2[N + 1][4].lower() and up(dz + 2 * eta0) < up(R)
        if ok:
            schritte[N] = it
            break
        F, J = F_J(X, gk1, N + 1)                                  # exakter Newton-Schritt, in Intervallen eingeschlossen
        d = J.solve(F)
        X = [X[k] - d[k, 0] for k in range(N + 2)]
    else:
        raise AssertionError(f"Übergang {N} -> {N + 1} nicht bestätigt")
print("   Intervall-Newton-Schritte vor Kantorovich je Übergang N->N+1:", schritte)
# N ≥ n0: Störungsargument um q_{N+1} mit monotonen Schranken (schlechtester Fall N = n0)
d = rmG + DREI ** (-(n0 + 2))
t = bG * gamma * d
assert up(t) < 1
b0 = bG / (1 - t)
eta0 = b0 * (resid(n0 + 1) + Lip_F * d)
h0, rm0, _ = kanto(b0, eta0)
assert up(d + 2 * eta0) < up(R) and up(d + rm0) < rpG.lower()
print(f"Teil 3: Newton vom aufgefüllten Vorgänger p_N* landet für jedes N ≥ 1 in p_(N+1)* (N < {n0} direkt, N ≥ {n0}: h0 ≤ {up(h0).str(4)}).")
print("ERGEBNIS (bewiesen, rechnergestützt): Der Zweig existiert für alle N, ‖p_N − T_N f‖_1 ≤ r_−(N) = O(6^−N), p_N → f auf |x| ≤ 1.")
