"""Gegenprobe: für eine andere reelle Tetration T(x) = sexp(x + ε sin 2πx) gelten die Nullmoden-Identitäten nicht.
Θ_T(x) = h(T x) − x, 1+Θ_T' = T'(x) h'(T x).  Aufruf: python s2b_kontrolle.py DPS"""
import sys, mpmath as mp
import s2_nullmoden as S
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = dps
L = S.L; PER = S.PER

def werte(eps, k):
    pts = S.knoten(k)
    th, lg = [], []
    for dl, dr, w in pts:
        x = dl if dl <= 0.5 else -dr
        u = x + eps * mp.sin(2 * mp.pi * x)
        s1, ds = S.sexp_reihe(u)
        l1 = mp.log1p(s1)
        hv, hd = S.B.h_mit_ableitung(l1)
        th.append(1 + hv - x)
        lg.append(mp.log(ds * (1 + 2 * mp.pi * eps * mp.cos(2 * mp.pi * x)) * hd / (1 + s1)))
    mid = len(pts) // 2
    th = S.entwirren(th[mid::-1], PER)[::-1][:-1] + S.entwirren(th[mid:], PER)
    lg = S.entwirren(lg[mid::-1], 2j * mp.pi)[::-1][:-1] + S.entwirren(lg[mid:], 2j * mp.pi)
    W = [p[2] for p in pts]
    M = lambda f: mp.fsum(w * f(i) for i, w in enumerate(W))
    mt = M(lambda i: th[i])
    return {"var": M(lambda i: th[i] ** 2) - mt ** 2,
            "eL": M(lambda i: mp.exp(L * th[i])) * M(lambda i: mp.exp(-L * th[i])) - 1,
            "log": M(lambda i: lg[i]),
            "inv": M(lambda i: mp.exp(-lg[i])) - 1}

for eps in (0, mp.mpf("1e-4"), mp.mpf("1e-3"), mp.mpf("1e-2")):
    r = werte(eps, 6)
    print(f"eps={mp.nstr(eps,3):6s}  " + "  ".join(f"{k}: {mp.nstr(abs(v),3):9s}" for k, v in r.items()), flush=True)
