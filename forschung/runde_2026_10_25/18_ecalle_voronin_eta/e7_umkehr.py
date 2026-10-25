# Koeffizienten der Umkehrung h^{-1} der oberen Hornabbildung, normierungsfreie Verhaeltnisse
import mpmath as mp, json
mp.mp.dps = 50
hi = json.load(open("cache/a_90_60_0.03_48_2.0.json"))
lo = json.load(open("cache/e1_40_32.json"))["1.5"]
a = {n: mp.mpc(*lo[n]) for n in range(1, 16)}
for n in range(1, 6): a[n] = mp.mpc(*hi[str(n)])
c = -1j * mp.pi / 3
def h(W): q = mp.expjpi(2 * W); return W + c + mp.fsum(a[n] * q**n for n in a)
Y = mp.mpf(3); M = 16
vals = []
for j in range(M):
    V = mp.mpc(mp.mpf(j) / M, Y); W = V - c
    for _ in range(60): W = W - (h(W) - V)  # Fixpunkt (h' ~ 1)
    vals.append(W - V + c)
inv = {}
for n in range(1, 5):
    inv[n] = mp.fsum(vals[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M * mp.exp(2 * mp.pi * n * Y)
print('Umkehrung: aI_1 =', mp.nstr(inv[1], 20))
for n in [2, 3, 4]:
    print(f' aI_{n}/aI_1^{n} =', mp.nstr(inv[n] / inv[1]**n, 15), f'   a_{n}/a_1^{n} =', mp.nstr(a[n] / a[1]**n, 12))
print(' 2 pi i - a2/a1^2 =', mp.nstr(2j * mp.pi - a[2] / a[1]**2, 15))
