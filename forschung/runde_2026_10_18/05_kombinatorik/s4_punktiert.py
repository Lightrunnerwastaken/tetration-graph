"""K5-Dreieck (Iterierte von x e^x) = Ketten im Verband der punktierten Partitionen (Chapoton–Vallette).

Punktierte Partition: Mengenpartition von [n], jeder Block hat ein ausgezeichnetes Element.
Ordnung: pi <= sigma  gdw. jeder sigma-Block Vereinigung von pi-Blöcken ist und sein Punkt einer der pi-Punkte ist.
Behauptung: n![x^n] (x e^x)^{o t} = sum_j B'(n,j) C(t,j),  B'(n,j) = #{0 < s_1 < ... < s_j, s_j maximal (ein Block)}.
Folgerungen: mu(0, M) = (-1)^{n-1} n^{n-2} für jedes Maximum M; Diagonale n!(n-1)!.
Brute force für n <= 6, gegen exakte Iterierte (Brüche) geprüft.
"""
import json
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, factorial
from pathlib import Path

HIER = Path(__file__).resolve().parent
NMAX = 6


def set_partitions(s):
    s = list(s)
    if not s:
        yield []
        return
    first, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [p[i] | {first}] + p[i + 1:]
        yield p + [{first}]


def pointed(n):
    out = []
    for p in set_partitions(range(n)):
        def rec(i, acc):
            if i == len(p):
                out.append(frozenset(acc))
                return
            for pt in p[i]:
                rec(i + 1, acc + [(frozenset(p[i]), pt)])
        rec(0, [])
    return out


def leq(a, b):
    if a == b:
        return True
    for blk, pt in b:
        sub = [(bb, pp) for bb, pp in a if bb <= blk]
        if frozenset().union(*[bb for bb, _ in sub]) != blk:
            return False
        if pt not in [pp for _, pp in sub]:
            return False
    return True


# exakte Iterierte von x e^x
def series_iter(N, k):
    f = [Fr(0)] + [Fr(1, factorial(m - 1)) for m in range(1, N + 1)]
    def mul(a, b):
        c = [Fr(0)] * (N + 1)
        for i, ai in enumerate(a):
            if ai:
                for j in range(N + 1 - i):
                    c[i + j] += ai * b[j]
        return c
    def comp(g):
        r = [Fr(0)] * (N + 1); p = [Fr(1)] + [Fr(0)] * N
        for m in range(N + 1):
            if f[m]:
                r = [ri + f[m] * pi for ri, pi in zip(r, p)]
            p = mul(p, g)
        return r
    g = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(k):
        g = comp(g)
    return g


res = {}
ok_all = True
for n in range(1, NMAX + 1):
    P = pointed(n)
    bottom = frozenset((frozenset([i]), i) for i in range(n))
    maxi = [p for p in P if len(p) == 1]
    up = {p: [q for q in P if q != p and leq(p, q)] for p in P}
    # Ketten von bottom aus zählen, die in einem Maximum enden, nach Länge
    from functools import lru_cache
    @lru_cache(None)
    def chains(p):  # dict Länge -> Anzahl Ketten p < s_1 < ... < s_j mit s_j maximal
        d = {}
        if len(p) == 1:
            d[0] = 1
        for q in up[p]:
            for j, c in chains(q).items():
                d[j + 1] = d.get(j + 1, 0) + c
        return d
    ch = chains(bottom)
    B = [ch.get(j, 0) for j in range(1, n)] if n > 1 else [1]
    if n == 1:
        B = [1]   # t-unabhängig: x
    # Vergleich mit Iterierten t = 0..n
    ok = True
    for t in range(0, n + 1):
        val = series_iter(n, t)[n] * factorial(n)
        pred = sum(B[j - 1] * comb(t, j) for j in range(1, len(B) + 1)) if n > 1 else 1
        ok &= (val == pred)
    # Möbius mu(0, M) über Philip Hall: sum_j (-1)^j #Ketten 0<...<M
    @lru_cache(None)
    def mu_to(p, M):
        if p == M:
            return 1
        return -sum(mu_to(q, M) for q in up[p] if leq(q, M))
    mus = sorted(set(mu_to(bottom, M) for M in maxi)) if n > 1 else [1]
    print(f"n={n}: #punktierte Partitionen={len(P)}, B'(n,j)={B}, Iterierte passen: {ok}, mu(0,M)={mus}, "
          f"(-1)^(n-1) n^(n-2)={(-1)**(n-1)*n**(n-2) if n>1 else 1}")
    ok_all &= ok and (n == 1 or mus == [(-1) ** (n - 1) * n ** (n - 2)])
    res[n] = {"anzahl": len(P), "B": B, "mu": mus}
print("alles bestätigt:", ok_all)
json.dump({"ergebnis": res, "ok": ok_all}, open(HIER / "daten" / "s4_punktiert.json", "w"), indent=1)
