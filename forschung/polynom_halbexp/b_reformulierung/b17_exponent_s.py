"""b17: Komplexer Exponent a(s) in c_{-1}(u_s) ~ N^{a(s)} für u_s = (1 − z/L)^s + (1 − z/L̄)^{s̄}.
Vermutung aus dem nichtlinearen Befund: a(s) = 2 − γ − s (generisch) oder nur am Kern-Exponenten s = 1+γ besonders?
Komplexer lokaler Exponent: log(c(N2)/c(N1))/log(N2/N1), gemittelt über eine Oszillationsperiode (N = 60 … 180)."""
import json, sys, math
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
ss = [1 + gam, mp.mpf(1) + mp.re(gam), 1 + gam + mp.mpf("0.5"), 1 + gam - mp.mpf("0.5"), mp.mpc(1 + mp.re(gam), 2 * mp.im(gam)), mp.mpf("3.0")]
tests = {}
for s in ss:
    out, cL, cLb = [], mp.mpc(1), mp.mpc(1)
    for j in range(601):
        out.append(mp.re(cL + cLb))
        cL = cL * (s - j) / (j + 1) * (-1 / L)
        cLb = cLb * (mp.conj(s) - j) / (j + 1) * (-1 / mp.conj(L))
    tests[mp.nstr(s, 4)] = (out, lambda z, s=s: mp.re((1 - z / L) ** s + (1 - z / mp.conj(L)) ** mp.conj(s)), s)
Ns = list(range(60, 181, 6))
res = {name: {} for name in tests}
for N in Ns:
    pw, g, Ji = baue(N, J=600)
    for name, (uj, uf, s) in tests.items():
        co = proj(N, pw, g, Ji, uj)
        res[name][N] = moden(N, co, uf)[-1]
for name, (uj, uf, s) in tests.items():
    d = res[name]
    # komplexer Fit: log c = a log N + b (kleinste Quadrate, komplex; Logarithmus stetig fortgesetzt)
    lc, prev = [], None
    for N in Ns:
        v = complex(mp.log(d[N]))
        if prev is not None:
            while v.imag - prev.imag > math.pi: v -= 2j * math.pi
            while v.imag - prev.imag < -math.pi: v += 2j * math.pi
        lc.append(v); prev = v
    xs = [math.log(N) for N in Ns]
    n = len(xs); mx = sum(xs) / n; my = sum(lc) / n
    a = sum((x - mx) * (y - my) for x, y in zip(xs, lc)) / sum((x - mx) ** 2 for x in xs)
    pred = complex(2 - gam - s)
    pred2 = complex(1 - gam - s)
    print(f"s={name}: a(s) Fit = {a.real:+.3f}{a.imag:+.3f}i   2−γ−s = {pred.real:+.3f}{pred.imag:+.3f}i   1−γ−s = {pred2.real:+.3f}{pred2.imag:+.3f}i")
