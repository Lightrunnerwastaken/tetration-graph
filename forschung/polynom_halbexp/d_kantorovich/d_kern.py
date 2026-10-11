"""Strang D — Kernfunktionen: Polynom-Wurzeln p_N hochgenau, F, J, Ränder (Bordering), Normen.

p_N(p_N(x)) = e^x + O(x^{N+1}).  J[k][j] = [x^k] p^j + [x^{k-j}] p'∘p.
Daten: gemeinsam/polynome_bis_194.json (60 Stellen) -> Newton-Verfeinerung in arb (approx-Löser) -> Zwischenspeicher hier.
"""
import json
import os
import pickle
import sys

from flint import arb, arb_mat, arb_poly, ctx

HIER = os.path.dirname(os.path.abspath(__file__))
GEMEINSAM = os.path.join(HIER, "..", "gemeinsam")
L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"


def setze_prec(bits):
    ctx.prec = bits


def potenzen(a, M):
    """[p^0, ..., p^M] modulo x^{M+1}  (p = Polynom mit Koeffizienten a)."""
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(M):
        P.append((P[-1] * f).truncate(M + 1))
    return P


def FJ(a, N, P=None):
    """F (Liste, Länge N+1) und J (arb_mat (N+1)x(N+1)) des Systems vom Grad N an der Stelle a (Länge N+1)."""
    if P is None:
        P = potenzen(a, N)
    g, fpf = arb_poly([0]), arb_poly([0])
    for m in range(N + 1):
        g += a[m] * P[m]
        if m:
            fpf += m * a[m] * P[m - 1]
    F, fak = [], arb(1)
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        fak = fak * k if k else fak
        F.append(g[k] - 1 / fak)
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def mitte_vec(v):
    return [arb(x.mid()) for x in v]


def loese(J, rhs):
    """approx. Lösung J x = rhs (rhs Liste) -> Liste."""
    n = len(rhs)
    B = arb_mat(n, 1, rhs)
    X = J.solve(B, algorithm="approx")
    return [X[i, 0] for i in range(n)]


def newton(a, N, tol_bits=None, maxit=30):
    tol = arb(2) ** (-(tol_bits or ctx.prec - 350))
    for _ in range(maxit):
        F, J = FJ(a, N)
        d = loese(J, F)
        a = [arb((a[i] - d[i]).mid()) for i in range(N + 1)]
        if max(abs(x.mid()) for x in d) < tol:
            return a
    raise RuntimeError(f"Newton ohne Konvergenz bei N={N}")


def lade_daten():
    d = json.load(open(os.path.join(GEMEINSAM, "polynome_bis_194.json")))
    return {int(k): v for k, v in d.items()}


def loesungen(nmax, bits=1200, cache=True):
    """{N: [arb a_0..a_N]} mit ~bits Genauigkeit; verfeinert aus den 60-Stellen-Daten; Zwischenspeicher als Strings."""
    datei = os.path.join(HIER, f"cache_p_{bits}.pkl")
    alt = {}
    if cache and os.path.exists(datei):
        alt = pickle.load(open(datei, "rb"))
    ctx.prec = bits
    daten = None
    aus, neu = {}, False
    digits = int(bits * 0.30103) + 5
    for N in range(1, nmax + 1):
        if N in alt:
            aus[N] = [arb(s) for s in alt[N]]
            continue
        if daten is None:
            daten = lade_daten()
        if N in daten:
            a0 = [arb(s) for s in daten[N]]
        else:
            a0 = aus[N - 1] + [arb(0)]
        aus[N] = newton(a0, N)
        alt[N] = [x.mid().str(digits, radius=False) for x in aus[N]]
        neu = True
        if cache and N % 20 == 0:
            pickle.dump(alt, open(datei, "wb"))
    if neu and cache:
        pickle.dump(alt, open(datei, "wb"))
    return aus


def L_wert():
    from flint import acb
    return acb(arb(L_RE), arb(L_IM))


def gew_norm(v, rho):
    """Σ |v_k| rho^k (arb)."""
    rho = arb(rho)
    s, w = arb(0), arb(1)
    for x in v:
        s += abs(x) * w
        w *= rho
    return s


def to_float(x):
    return float(x.mid())


if __name__ == "__main__":
    import time
    sys.stdout.reconfigure(encoding="utf-8")
    t = time.time()
    p = loesungen(int(sys.argv[1]) if len(sys.argv) > 1 else 40, bits=int(sys.argv[2]) if len(sys.argv) > 2 else 1200)
    print("Zeit", time.time() - t)
    N = max(p)
    print(N, p[N][0].str(40), p[N][N].str(20))
