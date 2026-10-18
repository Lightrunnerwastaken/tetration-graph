"""Gegenprobe zu E2 mit anderem Weg: T(x) = sexp_K(φ(x)) punktweise (Reihe um nächstes Zentrum bzw. log-Weg), Ableitung per mp.diff."""
import mpmath as mp, rk_lib as R, stoerungen as S
R.init(80)
def sexpK(u):
    if u < mp.mpf('-0.5'):
        return mp.log(sexpK(u + 1))
    c0 = min(R._Z, key=lambda c: abs(u - c)); co = R._Z[c0]
    return mp.polyval(co[::-1], u - c0)
for name, x, k in [("sin:0.001", '-0.6', 5), ("sin:0.03", '0', 3), ("sin:0.03", '-1.2', 3), ("kneser", '-0.6', 5)]:
    st = S.katalog(name)
    f = lambda t: sexpK(t + mp.re(S.p_komplex(st, t)[0]))
    d = mp.diff(f, mp.mpf(x), k, h=mp.mpf('1e-8'))
    mp.mp.dps = 50; c = R.tet(st, mp.mpf(x), k + 1)[k] * mp.factorial(k); mp.mp.dps = 80
    print(name, x, "T^(%d):" % k, mp.nstr(d, 12), " Reihen-Weg:", mp.nstr(c, 12))
