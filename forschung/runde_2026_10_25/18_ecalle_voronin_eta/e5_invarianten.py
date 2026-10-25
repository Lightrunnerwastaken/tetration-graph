import mpmath as mp, json
mp.mp.dps = 44
d = json.load(open("cache/a_90_60_0.03_48_2.0.json"))
a = {int(n): mp.mpc(*v) for n, v in d.items()}
a1 = a[1]
print('|a1| =', mp.nstr(abs(a1), 40)); print('arg a1 =', mp.nstr(mp.arg(a1), 40), ' /pi =', mp.nstr(mp.arg(a1) / mp.pi, 30))
print('|a1|^2 =', mp.nstr(abs(a1)**2, 40))
print('ln|a1| =', mp.nstr(mp.log(abs(a1)), 40))
for n in [2, 3, 4, 5]:
    r = a[n] / a1**n
    print(f'a{n}/a1^{n} =', mp.nstr(r, 30), ' |.| =', mp.nstr(abs(r), 20))
print('|a_n|^(1/n):', [mp.nstr(abs(a[n]) ** (mp.mpf(1) / n), 8) for n in range(1, 6)])
x = abs(a1)
for name, v in [('|a1|', x), ('|a1| e^{-pi^2/3}', x * mp.exp(-mp.pi**2 / 3)), ('|a1|/(2pi)', x / (2 * mp.pi)), ('|a1| e^{-pi^2/3}/(2pi)', x * mp.exp(-mp.pi**2/3)/(2*mp.pi))]:
    print(name, mp.nstr(v, 30), mp.identify(v, ['pi', 'e', 'log(2)', 'log(pi)', 'euler']), mp.identify(mp.log(v), ['pi', 'pi**2', 'log(2)', 'log(pi)', 'euler', 'log(3)']))
