import mpmath as mp, json
mp.mp.dps = 30
lo = json.load(open("cache/e1_40_32.json"))["1.5"]
a = [mp.mpc(*lo[n]) for n in range(16)]
for n in range(1, 15):
    print(n, mp.nstr(abs(a[n]) ** (mp.mpf(1) / n), 8), mp.nstr(abs(a[n + 1] / a[n]), 8), ' arg(a_{n+1}/a_n)/pi =', mp.nstr(mp.arg(a[n + 1] / a[n]) / mp.pi, 8))
eta = mp.e ** (1 / mp.e)
print('A = 2 pi^2 sqrt(eta/(2e)) =', mp.nstr(2 * mp.pi**2 * mp.sqrt(eta / (2 * mp.e)), 12))
print('e^{-2pi^2/3} =', mp.nstr(mp.exp(-2 * mp.pi**2 / 3), 15), '  e^{-pi^2/3} =', mp.nstr(mp.exp(-mp.pi**2 / 3), 15))
