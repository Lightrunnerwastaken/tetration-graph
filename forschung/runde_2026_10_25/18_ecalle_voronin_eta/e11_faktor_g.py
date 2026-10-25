# Universeller Faktor g(k) aus kl.verm.faktor gegen die Hornabbildung: g(k) =? |aI_k| / |aI_1|  (aI = Koeff. der Umkehrung h^{-1})
import mpmath as mp, json
mp.mp.dps = 40
hi = json.load(open("cache/a_90_60_0.03_48_2.0.json")); lo = json.load(open("cache/e1_40_32.json"))["1.5"]
a = {n: mp.mpc(*lo[n]) for n in range(1, 16)}
for n in range(1, 6): a[n] = mp.mpc(*hi[str(n)])
c = -1j * mp.pi / 3
def h(W): q = mp.expjpi(2 * W); return W + c + mp.fsum(a[n] * q**n for n in a)
Y = mp.mpf('2.2'); M = 32
vals = []
for j in range(M):
    V = mp.mpc(mp.mpf(j) / M, Y); W = V - c
    for _ in range(80): W = W - (h(W) - V)
    vals.append(W - V + c)
aI = {n: mp.fsum(vals[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M * mp.exp(2 * mp.pi * n * Y) for n in range(1, 11)}
d = json.load(open("../../../experimente/ergebnisse/exp.sqrt2_abel.json", encoding="utf-8"))["daten"]["basen"]
print("k  |aI_k|/|aI_1|      k^-5/4     " + "  ".join("r_k(b=" + B["basis"][:6] + ")" for B in d))
for k in range(1, 11):
    g = abs(aI[k]) / abs(aI[1])
    print(k, mp.nstr(g, 10), mp.nstr(mp.mpf(k) ** (-1.25), 6), "  ".join(str(B["r_k"][k - 1])[:10] if len(B["r_k"]) >= k else "-" for B in d))
print('|aI_1| =', mp.nstr(abs(aI[1]), 20), ' = |a1| e^{-2pi^2/3} =', mp.nstr(abs(a[1]) * mp.exp(-2 * mp.pi**2 / 3), 20))
json.dump({k: [str(mp.re(v)), str(mp.im(v))] for k, v in aI.items()}, open("cache/umkehr_aI.json", "w"))
