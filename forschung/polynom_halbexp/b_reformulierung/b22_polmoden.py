"""b22 (Runde 2): Modus-Antwort auf einen Pol. m_N(ζ) := Modus-(−1)-Koeffizient von θ_v, v = Π_N k_ζ − k_ζ, k_ζ(z) = 1/(ζ − z).
Exakt (Cauchy): c(u) = (1/2πi)∮ u(ζ) m_N(ζ) dζ über jede Kontur um 0, a_0 im Holomorphiegebiet von u.
Normierung nach der „+1-Regel“: K_N(ζ) := m_N(ζ) · [z^N](1−z/L)^{1−γ} / (N · ζ^{−N−1}).
Frage: Ist K_N(L(1 + v/N)) eine feste Funktion K(v) (innere Skala 1/N)? Dann ist K(s) (b19) ihre Mellin-Transformierte."""
import json, sys
import mpmath as mp
from flint import arb, acb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])


def rhs_pol(N, pw, g, zeta):
    """T_N[1/(ζ − f(x)) + g(x)/(ζ − x)] als komplexe Koeffizientenliste (mpmath)."""
    a = [mp.mpf(c) for c in K[:N + 1]]
    c0 = 1 / (zeta - a[0])
    s1 = [c0]
    for n in range(1, N + 1):
        s = mp.fsum(a[k] * s1[n - k] for k in range(1, n + 1))
        s1.append(c0 * s)
    gk = [mp.mpf(g[k].mid().str(100, radius=False)) for k in range(N + 1)]
    inv = [zeta ** (-(j + 1)) for j in range(N + 1)]
    s2 = [mp.fsum(gk[i] * inv[n - i] for i in range(n + 1)) for n in range(N + 1)]
    return [s1[n] + s2[n] for n in range(N + 1)]


def loese(Ji, b):
    br = arb_mat([[arb(mp.nstr(mp.re(x), 100))] for x in b]); bi = arb_mat([[arb(mp.nstr(mp.im(x), 100))] for x in b])
    xr = Ji * br; xi = Ji * bi
    n = len(b)
    return [mp.mpc(mp.mpf(xr[k, 0].mid().str(100, radius=False)), mp.mpf(xi[k, 0].mid().str(100, radius=False))) for k in range(n)]


def mode(coef, fun, n=-1):
    th = [(horner(coef, z) - fun(z)) * gw for (w, z), gw in zip(pts, gewicht)]
    return mp.fsum(t * mp.exp(-2j * mp.pi * n * w) for t, (w, z) in zip(th, pts)) / M


mp.mp.dps = 100
res = {}
vs = [0.25, 0.5, 1, 2, 4, 8, 16]
for N in (60, 100, 140):
    pw, g, Ji = baue(N, J=N)
    deltaN = mp.binomial(1 - gam, N) * (-1 / L) ** N
    zeile = []
    for v in vs:
        zeta = L * (1 + mp.mpf(v) / N)
        P = loese(Ji, rhs_pol(N, pw, g, zeta))
        m = mode(P, lambda z: 1 / (zeta - z))
        Kv = m * deltaN / (N * zeta ** (-N - 1))
        res.setdefault(N, {})[v] = Kv
        zeile.append(f"v={v}: {mp.nstr(Kv, 4)}")
    print(f"N={N}: K_N(L(1+v/N)) = " + " | ".join(zeile), flush=True)
json.dump({N: {str(v): [str(mp.re(x)), str(mp.im(x))] for v, x in d.items()} for N, d in res.items()}, open("b22_polmoden.json", "w"))
