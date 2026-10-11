"""b21: Ist der Zusatzfaktor N in der Auswahlregel strukturell? Deformiertes lineares Problem
    J'_λ δ := δ∘f + λ (f'∘f) δ   (f = f_K),   Π'_N u := (J'_N)⁻¹ T_N J'_λ u.
Kern: δ_θ mit θ(w + ½) = −λ θ(w)  ⇒  θ = λ^{2w}·(antiperiodisch), Konkurrenzmodus bei L: Exponent s_c = 1 − γ − 2 log λ / L.
Vorhersage (Regel „N × Koeffizientenquotient“): Amplitude ~ N^{1 + s_c − s};  ohne Zusatzfaktor: N^{s_c − s}."""
import json, sys, math
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("def baue")[0])


def baue_l(N, lam, J=600):
    f = arb_poly([arb(c) for c in K[:N + 1]])
    pw = [arb_poly([1])]
    for _ in range(J):
        pw.append((pw[-1] * f).truncate(N + 1))
    ex_ser = [arb(1)]
    for k in range(1, N + 1):
        ex_ser.append(ex_ser[-1] / k)
    fp = [arb(K[k]) * k for k in range(1, N + 2)]
    gk = []
    for n in range(N + 1):
        s = ex_ser[n] - sum((fp[k] * gk[n - k] for k in range(1, n + 1)), arb(0))
        gk.append(s / fp[0])
    g = arb_poly([arb(lam) * x for x in gk])          # λ·(f'∘f) = λ e^x / f'(x)
    Jm = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Jm[k, j] = pw[j][k] + (g[k - j] if j <= k else 0)
    return pw, g, Jm.inv()


def proj_l(N, pw, g, Ji, uj):
    rhs = [arb(0)] * (N + 1)
    for j in range(len(pw)):
        cj = arb(mp.nstr(uj[j], 70))
        for k in range(N + 1):
            rhs[k] += cj * pw[j][k]
    gu = (g * arb_poly([arb(mp.nstr(uj[j], 70)) for j in range(N + 1)])).truncate(N + 1)
    x = Ji * arb_mat([[rhs[k] + gu[k]] for k in range(N + 1)])
    return [mp.mpf(x[k, 0].mid().str(70, radius=False)) for k in range(N + 1)]


def mode_m1(coef, ufun, lam):
    th = [(horner(coef, z) - ufun(z)) * gw * mp.power(lam, -2 * w) for (w, z), gw in zip(pts, gewicht)]
    return mp.fsum(t * mp.exp(2j * mp.pi * w) for t, (w, z) in zip(th, pts)) / M


s = mp.mpf("5.447")
uj = []
cL, cLb = mp.mpc(1), mp.mpc(1)
for j in range(601):
    uj.append(mp.re(cL + cLb))
    cL = cL * (s - j) / (j + 1) * (-1 / L)
    cLb = cLb * (s - j) / (j + 1) * (-1 / mp.conj(L))
ufun = lambda z: mp.re((1 - z / L) ** s + (1 - z / mp.conj(L)) ** s)
Ns = list(range(60, 181, 12))
for lam in (mp.mpf(1), mp.mpf(2), mp.mpf("0.5")):
    vals = {}
    for N in Ns:
        pw, g, Ji = baue_l(N, lam)
        vals[N] = mode_m1(proj_l(N, pw, g, Ji, uj), ufun, lam)
    lc, prev = [], None
    for N in Ns:
        v = complex(mp.log(vals[N]))
        if prev is not None:
            while v.imag - prev.imag > math.pi: v -= 2j * math.pi
            while v.imag - prev.imag < -math.pi: v += 2j * math.pi
        lc.append(v); prev = v
    xs = [math.log(N) for N in Ns]
    n = len(xs); mx = sum(xs) / n; my = sum(lc) / n
    a = sum((x - mx) * (y - my) for x, y in zip(xs, lc)) / sum((x - mx) ** 2 for x in xs)
    sc = 1 - gam - 2 * mp.log(lam) / L
    p1 = complex(1 + sc - s); p0 = complex(sc - s)
    print(f"λ = {mp.nstr(lam,3)}: Fit a = {a.real:+.3f}{a.imag:+.3f}i;  1+s_c−s = {p1.real:+.3f}{p1.imag:+.3f}i;  s_c−s = {p0.real:+.3f}{p0.imag:+.3f}i;"
          f"  |c| bei N=60/180: {mp.nstr(abs(vals[60]),3)}/{mp.nstr(abs(vals[Ns[-1]]),3)}", flush=True)
