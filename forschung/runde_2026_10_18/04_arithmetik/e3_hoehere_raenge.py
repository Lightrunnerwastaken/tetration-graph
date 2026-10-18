"""E3: hoehere Raenge.
 (a) Bewertungsformel bei grossen Hoehen H (100, 400, 1000)
 (b) 'Frischer Block': F_j(H) = ((a^^H - alpha)/p^{v_H}) mod p^j haengt fuer grosse H nur von H mod ord(u, p^j) ab
     (u = Einheitsteil von lambda = log_p alpha) -> bei Pentation a^^^n (Hoehe M = a^^^(n-1)) schliesslich konstant in n.
 (c) Tabelle stabiler Dezimalstellen fuer Pentation/Hexation (Formel)
"""
from ar_lib import *

print("(a) grosse Hoehen")
for p, a in ((5, 2), (5, 3), (3, 2), (7, 3), (2, 3), (2, 5), (13, 6)):
    for H in (100, 300):
        if p == 2:
            pred = vp(a - 1, 2) + H * (vp(a * a - 1, 2) - 1)
        else:
            V = vp(pow(a, p - 1) - 1, p); t = stab_height(a, order_mod(a, p)); pred = V * (H - t)
        K = pred + 5
        q = p ** K
        seq = tower_run(a, p, K, 8, H + K + 10)
        al = seq[-1]
        obs = vp((seq[H - 8] - al) % q, p)
        print(f"  p={p} a={a} H={H}: beobachtet {obs}, Formel {pred}", "ok" if obs == pred else "FEHLER")

print("(b) frischer Block")
for p, a, j in ((5, 2, 3), (5, 3, 3), (3, 2, 4), (7, 3, 2), (2, 3, 5), (2, 5, 5)):
    N = 150
    q = p ** N
    seq = tower_run(a, p, N, 6, 6 + N + 10)
    al = seq[-1]
    lam = al * padic_log(a, p, N) % q
    V = vp(lam, p); u = (lam // p ** V) % p ** j
    o = 1; x = u
    while x % p ** j != 1:
        x = x * u % p ** j; o += 1
    blocks = {}
    for H in range(6, 200):
        d = (seq[H - 6] - al) % q
        v = vp(d, p)
        if v + j > N:
            break
        blocks.setdefault(H % o, set()).add((d // p ** v) % p ** j)
    # Ab welcher Hoehe gilt die Abhaengigkeit nur von H mod o? Pruefe alle H >= 12
    blocks2 = {}
    for H in range(12, 200):
        d = (seq[H - 6] - al) % q
        v = vp(d, p)
        if v + j > N:
            break
        blocks2.setdefault(H % o, set()).add((d // p ** v) % p ** j)
    print(f"  p={p} a={a} j={j}: ord(u mod p^j)={o}; Bloecke je Restklasse H mod {o} (H>=12):",
          {k: sorted(vv) for k, vv in sorted(blocks2.items())} if all(len(vv) == 1 for vv in blocks2.values()) else "NICHT eindeutig")

print("(c) stabile Dezimalstellen D(a^^M) = min(v_2, v_5)")
def vform(a, p, M):
    if a % p == 0:
        e = tower_exact(a, M - 1)
        return float('inf') if e is None else vp(a, p) * e
    if p == 2:
        return vp(a - 1, 2) + M * (vp(a * a - 1, 2) - 1)
    V = vp(pow(a, p - 1) - 1, p); t = stab_height(a, order_mod(a, p))
    return V * (M - t) if M >= t else None
# direkte Kontrolle der Dezimalformel fuer M <= 40
bad = 0
for a in range(2, 60):
    if a % 10 == 0:
        continue
    al = tower_mod(a, 400, 10 ** 60)
    for M in range(3, 40):
        x = tower_mod(a, M, 10 ** 60)
        d = 0
        while d < 60 and (x - al) % 10 ** (d + 1) == 0:
            d += 1
        f = min(vform(a, 2, M), vform(a, 5, M))
        if f is None or f >= 60:
            continue
        if d != f:
            bad += 1
print("  Dezimalformel M in [3,40), a in [2,60): Abweichungen", bad)
for a in (2, 3, 4, 5, 6, 7, 11, 12, 13):
    s = []
    for M, name in ((a, f"{a}^^^2 (M={a})"),):
        pass
    M2 = a  # a^^^2 = a^^a
    M3 = tower_exact(a, a)  # a^^^3 = a^^(a^^a)
    f2 = min(vform(a, 2, M2), vform(a, 5, M2))
    f3 = None if M3 is None else min(vform(a, 2, M3), vform(a, 5, M3))
    V10 = min((vp(a * a - 1, 2) - 1) if a % 2 else float('inf'), vp(pow(a, 4) - 1, 5) if a % 5 else float('inf'))
    print(f"  a={a}: V10={V10}; D(a^^^2)={f2}; D(a^^^3)={f3 if f3 is not None else 'M zu gross: ' + 'V10*(a^^a) + c'}")
