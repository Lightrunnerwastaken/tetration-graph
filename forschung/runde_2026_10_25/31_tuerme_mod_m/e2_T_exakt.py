"""E2: exakte maximale Stabilisierungshoehe T(m) = max_a t(a,m).

(i)  Brute force ueber ALLE a (ein volles Restsystem mod Lam(m) oberhalb von A0, dazu alle a < A0) fuer kleine m.
(ii) Formel aus Satz R: T(m) = max( l(m)+delta(m), max_{a<A0} t(a,m) ),
     delta(q) fuer Primzahlpotenzen: delta(2^e)=0 (e<=3), -1 (e>=4); delta(p^k)=delta(p) (p ungerade);
     p ungerade gut (delta=0)  <=>  eine Komponente s^e || p-1 mit l(s^e)=l(p-1) ist gut.
     delta(m) = 0 <=> eine Komponente q || m mit l(q) = l(m) ist gut, sonst -1.
"""
import sys, time
from functools import lru_cache
from collections import Counter
from tlib import fac, ell, Lam, t_rec

A0 = 30


@lru_cache(maxsize=None)
def delta_pp(p, k):
    if p == 2:
        return 0 if k <= 3 else -1
    return delta_m(p - 1) if False else _delta_odd(p)


@lru_cache(maxsize=None)
def _delta_odd(p):
    return delta_m(p - 1)


def delta_m(m):
    if m == 1:
        return 0
    comps = fac(m)
    L = max(ell(p ** k) for p, k in comps)
    return 0 if any(ell(p ** k) == L and delta_pp(p, k) == 0 for p, k in comps) else -1


def T_formula(m):
    v = ell(m) + delta_m(m)
    for a in range(1, A0):
        v = max(v, t_rec(a, m))
    return v


def T_brute(m):
    L = Lam(m)
    best = -9
    for a in range(1, A0 + L):
        best = max(best, t_rec(a, m))
    return best


if __name__ == "__main__":
    MB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
    t0 = time.time(); ok = 0; skipped = 0
    for m in range(2, MB + 1):
        if Lam(m) > LMAX:
            skipped += 1; continue
        tb, tf = T_brute(m), T_formula(m)
        if tb != tf:
            print("ABWEICHUNG m =", m, tb, tf)
        else:
            ok += 1
    print(f"Brute force gegen Formel: {ok} Moduln gleich (m <= {MB}, Lam(m) <= {LMAX}), {skipped} uebersprungen, {time.time()-t0:.1f}s")
