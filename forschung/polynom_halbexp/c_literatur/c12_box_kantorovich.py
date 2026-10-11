"""Strang C, Runde 2, c12: Kantorovich in der Kastennorm |v|_m = max_k |v_k|/m_k (gewichtetes l^∞ mit Profil m).

Zentrum q = q_N (Kneser), A = J_N(q), δ1 = -A^{-1}F_N(q) (erster Newton-Schritt), G_{k,ij} = (A^{-1}F''(q)[e_i,e_j])_k.
  α_m = |δ1|_m,  ω_m = max_k Σ_{i,j} |G_{k,ij}| m_i m_j / m_k   (Supremum über den ganzen Kasten, mit Beträgen),
  h_m = α_m ω_m   (Kantorovich-Bedingung h ≤ 1/2, F'' an der Stelle q; Variation über den Kasten ist dritter Ordnung).
Profile:  P1: m = lokales Maximum von |δ1| (Fenster ±2) mit Boden 1e-8·max_l(|δ1_l|ρ^l)·ρ^{-k};
          P2: m_k = ρ^{-k}(1+k)^{-β}·c mit β aus einem Fit an |δ1| (Hüllkurve);
          P3: m_k = |f_k| + |δ1_k| (Kneser-Koeffizienten als Boden).
Ausgabe je N und Profil: α_m, ω_m, h_m, und das k mit dem grössten Beitrag.
"""
import json, sys, os, math
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf
from flint import arb_mat, arb, ctx

mp.dps = 120
ctx.prec = 500
HIER = os.path.dirname(os.path.abspath(__file__))
G0 = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G0, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
RHO = mpf("1.3745570107437075")


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def kompon(c, pot, n):
    r = [mpf(0)] * (n + 1)
    for j, cj in enumerate(c[: n + 1]):
        if cj == 0:
            continue
        for k in range(n + 1):
            r[k] += cj * pot[j][k]
    return r


def shift(v, s, n):
    return ([mpf(0)] * s + v)[: n + 1] if s <= n else [mpf(0)] * (n + 1)


def a2(x):
    return arb(mp.nstr(x, 110))


NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [20, 40, 60, 80, 100]
for N in NS:
    q = KN[: N + 1]
    n = N
    pot = [[mpf(1)] + [mpf(0)] * n]
    for _ in range(n):
        pot.append(mul(pot[-1], q, n))
    dq = [k * q[k] for k in range(1, n + 1)]
    ddq = [k * (k - 1) * q[k] for k in range(2, n + 1)]
    dpp = kompon(dq, pot, n)
    ddpp = kompon(ddq, pot, n)
    A = arb_mat([[a2(pot[j][k] + (dpp[k - j] if j <= k else 0)) for j in range(n + 1)] for k in range(n + 1)])
    Ai = A.inv()
    qq = kompon(q, pot, n)
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(n + 1)]
    d1 = Ai * arb_mat([[a2(-x)] for x in F])
    d1 = [float(abs(mpf(d1[k, 0].mid().str(100, radius=False)))) for k in range(n + 1)]
    paare, spalten = [], []
    for i in range(n + 1):
        for j in range(i, n + 1):
            v = [mpf(0)] * (n + 1)
            if i >= 1:
                t = shift(pot[i - 1], j, n)
                v = [x + i * y for x, y in zip(v, t)]
            if j >= 1:
                t = shift(pot[j - 1], i, n)
                v = [x + j * y for x, y in zip(v, t)]
            t = shift(ddpp, i + j, n)
            v = [x + y for x, y in zip(v, t)]
            paare.append((i, j))
            spalten.append(v)
    Bm = arb_mat([[a2(spalten[c][k]) for c in range(len(spalten))] for k in range(n + 1)])
    Gm = Ai * Bm
    Gabs = [[abs(float(Gm[k, c].mid())) for c in range(len(paare))] for k in range(n + 1)]
    rho = float(RHO)
    skal = max(d1[l] * rho ** l for l in range(n + 1))
    P1 = []
    for k in range(n + 1):
        lok = max(d1[max(0, k - 2): k + 3])
        P1.append(max(lok, 1e-8 * skal * rho ** (-k)))
    # P2: Hüllkurve ρ^{-k}(1+k)^{-β}: Fit von log(|δ1_k|ρ^k) gegen log(1+k) für k in [N/2, N]
    xs = [math.log(1 + k) for k in range(n // 2, n + 1)]
    ys = [math.log(max(d1[k] * rho ** k, 1e-300)) for k in range(n // 2, n + 1)]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    beta = -sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    c2 = max(d1[k] * rho ** k * (1 + k) ** beta for k in range(n + 1))
    P2 = [c2 * rho ** (-k) * (1 + k) ** (-beta) for k in range(n + 1)]
    P3 = [abs(float(KN[k])) + d1[k] for k in range(n + 1)]
    for name, m in (("P1", P1), ("P2", P2), ("P3", P3)):
        alpha = max(d1[k] / m[k] for k in range(n + 1))
        best, kb = 0.0, -1
        for k in range(n + 1):
            s = 0.0
            row = Gabs[k]
            for c, (i, j) in enumerate(paare):
                s += row[c] * m[i] * m[j] * (1 if i == j else 2)
            s /= m[k]
            if s > best:
                best, kb = s, k
        print(f"N={N:3d} {name}{'(beta=%.2f)' % beta if name == 'P2' else ''}: alpha_m={alpha:.3e} omega_m={best:.3e} "
              f"h_m={alpha * best:.3e} (max bei k={kb})")
    sys.stdout.flush()
