"""E1: Kneser — Vorzeichen aller Taylor-Koeffizienten c_k(x) (k < N) auf einem Gitter.
Aufruf: python e1_kneser_ungerade.py N dps schritt"""
import sys, json, mpmath as mp, rk_lib as R
N = int(sys.argv[1]); dps = int(sys.argv[2]); h = mp.mpf(sys.argv[3])
R.init(dps)
xs = []; x = mp.mpf('-1.95')
while x <= 2.0001: xs.append(x); x += h
res = []
min_odd = None
for x in xs:
    c = R.kneser(x, N)
    odd_neg = [k for k in range(1, N, 2) if c[k] <= 0]
    even_pos = [k for k in range(2, N, 2) if c[k] > 0]
    # normierte Grösse (-1)^{k-1} c_k k (x+2)^k  (→ 1 für k → ∞, wenn −2 dominiert)
    q = [(-1) ** (k - 1) * c[k] * k * (x + 2) ** k for k in range(1, N)]
    mo = min((c[k] * k * (x + 2) ** k, k) for k in range(1, N, 2))
    res.append({"x": mp.nstr(x, 6), "odd_neg": odd_neg, "even_pos_max": max(even_pos) if even_pos else 0,
                "min_odd_norm": [mp.nstr(mo[0], 6), mo[1]], "q_last": mp.nstr(q[-1], 8)})
    print(mp.nstr(x, 4), "ungerade<=0:", odd_neg, " gerade>0 bis k =", max(even_pos) if even_pos else 0,
          " min ungerade normiert", mp.nstr(mo[0], 5), "bei k", mo[1], " q_N", mp.nstr(q[-1], 6), flush=True)
json.dump(res, open(f"daten/e1_N{N}_dps{dps}.json", "w"), indent=0)
