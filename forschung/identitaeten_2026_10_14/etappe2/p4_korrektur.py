"""Zweite Ordnung des Modells: ℓ2·ln b = c + log u + δ1·u + O(u²),  u = 1 − q,
δ1 = 1/2 + iκ/(2π),  κ = B/A − A/2,  sexp(z) = 1 + A z + B z² + …
Korrigiertes Modell  Φ2 = Φ + δ1·u·Φ_λ  mit  u·dΦ/dλ = −u²·Φ'(q)  ⇒
b2_n = b_n − δ1·[(n+1)b_{n+1} − 2n b_n + (n−1) b_{n−1}].
Erwartung: |a_n − b2_n|/|a_n| = O(n^{-2}·langsam), statt O(n^{-1}) für b_n.
Aufruf: python p4_korrektur.py <basis e|2> <datei a_n> <datei modell>
"""
import json
import sys
from pathlib import Path

import mpmath as mp
from kern import sexp

B_, FA, FM = sys.argv[1], sys.argv[2], sys.argv[3]
mp.mp.dps = 60
basis = "exp(1)" if B_ == "e" else B_
K, rho = 64, mp.mpf("0.1")
pts = [rho * mp.expjpi(2 * mp.mpf(j) / K) for j in range(K)]
vals = sexp(basis, pts, 60)
co = [mp.fsum(v * mp.expjpi(-2 * mp.mpf(j * k) / K) for j, v in enumerate(vals)) / K / rho ** k for k in range(3)]
A, Bc = mp.re(co[1]), mp.re(co[2])
kappa = Bc / A - A / 2
d1 = mp.mpf(1) / 2 + 1j * kappa / (2 * mp.pi)
print("A =", mp.nstr(A, 25), " B =", mp.nstr(Bc, 25), " kappa =", mp.nstr(kappa, 20), " delta1 =", mp.nstr(d1, 20))
c = lambda v: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))
a = [c(v) for v in json.loads(Path(FA).read_text())["a"]]
b = [c(v) for v in json.loads(Path(FM).read_text())["b"]]
N = min(len(a), len(b)) - 2
b2 = [None] + [b[n] - d1 * ((n + 1) * b[n + 1] - 2 * n * b[n] + (n - 1) * b[n - 1]) for n in range(1, N + 1)]
zeilen = []
print("    n   |a-b|/|a|      |a-b2|/|a|     n^2 |a-b2|/|a|")
for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4000]:
    r1 = abs(a[n] - b[n]) / abs(a[n])
    r2 = abs(a[n] - b2[n]) / abs(a[n])
    zeilen.append([n, mp.nstr(r1, 5), mp.nstr(r2, 5), mp.nstr(n * n * r2, 5)])
    print(f"{n:5d}  {mp.nstr(r1,5):>12}  {mp.nstr(r2,5):>12}  {mp.nstr(n*n*r2,5):>12}")
Path("daten", f"korrektur_{B_}.json").write_text(json.dumps({"A": mp.nstr(A, 30), "B": mp.nstr(Bc, 30),
                                                            "delta1": mp.nstr(d1, 30), "zeilen": zeilen}, indent=1))
