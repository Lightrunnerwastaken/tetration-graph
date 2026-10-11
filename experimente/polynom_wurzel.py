"""Polynom-Wurzel von exp: p_N mit p_N(p_N(x)) = e^x + O(x^{N+1}), Zweig stetig von p_1 = x + 1/2 aus (Newton in arb, python-flint)."""
import mpmath as mp
from flint import arb, arb_mat, arb_poly, ctx


def _schritt(a, N):
    for _ in range(80):
        f = arb_poly(a)
        P = [arb_poly([1])]
        for _m in range(N):
            P.append((P[-1] * f).truncate(N + 1))
        g, fpf = arb_poly([0]), arb_poly([0])
        for m in range(N + 1):
            g += a[m] * P[m]
            if m:
                fpf += m * a[m] * P[m - 1]
        F, J, fak = arb_mat(N + 1, 1), arb_mat(N + 1, N + 1), arb(1)
        for k in range(N + 1):
            fak = fak * k if k else fak
            F[k, 0] = g[k] - 1 / fak
            for j in range(N + 1):                       # ∂/∂a_j von p(p(x)): p(x)^j + p'(p(x))·x^j
                J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
        d = J.solve(F)
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
        if max(abs(d[i, 0].mid()) for i in range(N + 1)) < arb(2) ** (-ctx.prec + 60):
            return a
    raise RuntimeError(f"Newton ohne Konvergenz bei N = {N}")


def loesungen(nmax, bits=500):
    """{N: [a_0, …, a_N]} als mpmath-Zahlen für N = 1 … nmax."""
    alt, ctx.prec = ctx.prec, bits
    try:
        aus, a = {}, [arb("0.5"), arb(1)]
        for N in range(1, nmax + 1):
            a = _schritt(a + [arb(0)] * (N + 1 - len(a)), N)
            aus[N] = [mp.mpf(c.mid().str(130, radius=False)) for c in a]
        return aus
    finally:
        ctx.prec = alt


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s
