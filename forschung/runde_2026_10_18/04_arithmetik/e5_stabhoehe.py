"""E5: Stabilisierungshoehe t(a,m) (kleinstes t>=-1 mit a^^j ≡ a^^t mod m fuer j>=t)
gegen die Laenge l(m) der Carmichael-Kette m -> lambda(m) -> ... -> 1."""
from ar_lib import stab_height, carmichael
from collections import Counter
def chainlen(m):
    l = 0
    while m > 1:
        m = carmichael(m); l += 1
    return l
res = Counter(); worst = []
for m in range(2, 2001):
    l = chainlen(m)
    T = max(stab_height(a, m) for a in range(2, 2 * m + 60))
    res[(T - l)] += 1
    if T > l: worst.append((m, T, l))
print("Verteilung T(m) - l(m) fuer m in [2,2000], a in [2, 2m+60):", sorted(res.items()))
print("Faelle T > l:", worst[:20])
