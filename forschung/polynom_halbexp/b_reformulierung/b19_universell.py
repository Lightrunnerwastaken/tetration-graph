"""b19: Ist die Modusamplitude universell durch den N-ten Taylorkoeffizienten bestimmt?
Test: R(s,N) := c_{-1}(u_s)·Γ(−s)·N^{s+γ−2}·(−L)^{-s}... d.h. c ≈ κ' N^{2−γ} L^N [z^N](1−z/L)^s ?
[z^N](1−z/L)^s = (−1/L)^N C(s,N) ≈ L^{−N} N^{−s−1}/Γ(−s). Wir vergleichen Q(s,N) := c_{-1}(u_s) / (N^{2−γ} L^N [z^N](1−z/L)^s)."""
import json, sys, math
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
ss = [1 + gam, mp.mpf(1) + mp.re(gam), 1 + gam + mp.mpf("0.5"), 1 + gam - mp.mpf("0.5"), mp.mpc(1 + mp.re(gam), 2 * mp.im(gam)), mp.mpf("2.5"), mp.mpf("3.5")]
tests = {}
for s in ss:
    out, cL, cLb = [], mp.mpc(1), mp.mpc(1)
    for j in range(601):
        out.append(mp.re(cL + cLb))
        cL = cL * (s - j) / (j + 1) * (-1 / L)
        cLb = cLb * (mp.conj(s) - j) / (j + 1) * (-1 / mp.conj(L))
    tests[mp.nstr(s, 4)] = (out, lambda z, s=s: mp.re((1 - z / L) ** s + (1 - z / mp.conj(L)) ** mp.conj(s)), s)
for N in (60, 90, 120, 150, 180):
    pw, g, Ji = baue(N, J=600)
    zeile = []
    for name, (uj, uf, s) in tests.items():
        co = proj(N, pw, g, Ji, uj)
        c = moden(N, co, uf)[-1]
        coefN = mp.binomial(s, N) * (-1 / L) ** N
        Q = c / (mp.power(N, 2 - gam) * L ** N * coefN)
        zeile.append(f"{name}: {mp.nstr(Q, 4)}")
    print(f"N={N}: " + " | ".join(zeile), flush=True)
