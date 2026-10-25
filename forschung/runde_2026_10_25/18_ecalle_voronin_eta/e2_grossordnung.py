# Grossordnung der formalen Fatou-Koeffizienten c_k gegen die Vorhersage aus a_1 (Resurgenz).
# Vorhersage: c_k ~ 2 Re[ D * i^{-k} * Gamma(k - 2 pi i/3) * (4 pi)^{-k + 2 pi i/3} ],  |D| = |a_1| e^{-pi^2/3} / (2 pi)
import mpmath as mp, sys, json
mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 120
K = int(sys.argv[1]) if len(sys.argv) > 1 else 120
N = K + 3
fac = [mp.mpf(1)]
for j in range(1, N + 3): fac.append(fac[-1] * j)
q = [1 / fac[j + 1] for j in range(N + 1)]
def mul(a, b):
    r = [mp.mpf(0)] * (N + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(N + 1 - i): r[i + j] += x * b[j]
    return r
def inv(a):
    r = [mp.mpf(0)] * (N + 1); r[0] = 1 / a[0]
    for n in range(1, N + 1):
        r[n] = -mp.fsum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r
iq = inv(q)
da = [q[k + 1] * (k + 1) for k in range(N)] + [mp.mpf(0)]
t = mul(da, iq)
lq = [mp.mpf(0)] + [t[k] / (k + 1) for k in range(N)]
total = [mp.mpf(0)] * (N + 1)
tt = [(-iq[k] if k else 1 - iq[0]) for k in range(N + 1)]
for k in range(N): total[k] += 2 * tt[k + 1]
for k in range(N + 1): total[k] += lq[k] / 3
total[0] -= 1
gser = [mp.mpf(0)] + q[:N]
c = [mp.mpf(0)] * (K + 1)
gk = gser[:]
for k in range(1, K + 1):
    ck = -total[k + 1] / (mp.mpf(k) / 2)
    c[k] = ck
    for j in range(N + 1): total[j] += ck * gk[j]
    total[k] -= ck
    gk = mul(gk, gser)
json.dump([str(x) for x in c[1:]], open(f"cache/c_{K}_{mp.mp.dps}.json", "w"))
a1 = mp.mpc('10.7488070931472942', '63.2462981521074514')
for k in list(range(10, K + 1, 10)) + [K - 1]:
    pred = lambda D: 2 * mp.re(D * mp.power(1j, -k) * mp.gamma(k - 2j * mp.pi / 3) * mp.power(4 * mp.pi, -k + 2j * mp.pi / 3))
    r = c[k] * (4 * mp.pi) ** k / mp.gamma(k)
    print(k, mp.nstr(c[k], 12), ' c_k(4pi)^k/Gamma(k) =', mp.nstr(r, 10))
