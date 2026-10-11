"""b11b: Vergleich linear (Π_N an f_K auf u_{1+γ}) mit nichtlinear (p_N, b10) und Exponenten-Fit über viele N."""
import json, sys, math
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
nl = {int(k): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in json.load(open("b10_amplitude.json")).items()}
tests = {}
for s in (mp.mpf("2.5"), 1 + gam, mp.mpf("5.5")):
    out, cL, cLb = [], mp.mpc(1), mp.mpc(1)
    for j in range(901):
        out.append(mp.re(cL + cLb))
        cL = cL * (s - j) / (j + 1) * (-1 / L)
        cLb = cLb * (mp.conj(s) - j) / (j + 1) * (-1 / mp.conj(L))
    tests[mp.nstr(s, 5)] = (out, lambda z, s=s: mp.re((1 - z / L) ** s + (1 - z / mp.conj(L)) ** mp.conj(s)))
res = {}
Ns = list(range(30, 191, 8))
for N in Ns:
    pw, g, Ji = baue(N, J=600)
    for name, (uj, uf) in tests.items():
        co = proj(N, pw, g, Ji, uj)
        res.setdefault(name, {})[N] = moden(N, co, uf)[-1]
    nlc = mp.conj(nl[N]) if N in nl else None
    lin = res[mp.nstr(1 + gam, 5)][N]
    print(f"N={N}: linear c(-1)[u_(1+γ)] = {mp.nstr(lin, 5)}, nichtlinear conj(c_N) = {mp.nstr(nlc, 5) if nlc is not None else '-'},"
          f" Quotient = {mp.nstr(nlc / lin, 5) if nlc is not None else '-'}", flush=True)
for name, d in res.items():
    xs = [math.log(N) for N in Ns]; ys = [math.log(float(abs(d[N]))) for N in Ns]
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    print(f"s={name}: Fit-Exponent von |c(-1)| über N=30..190: {sl:.3f}")
