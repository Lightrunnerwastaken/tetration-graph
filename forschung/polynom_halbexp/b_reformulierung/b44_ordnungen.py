"""b44 (Runde 4): Ordnungen der Entwicklung p_N = x₁ + d^{(2)} + d^{(3)} + … an echten Daten.
d^{(2)} = −½ J⁻¹T_N D²Φ(f)[u₁,u₁] = Π_N[u₁u₁'/(2f')] − (N+1)ρ^f u₁(0)/(2f'(0))·t_N   (Korollar B-Ξ₂), u₁ = x₁ − f_K.
Gemessen: ε_N = sup_D̄|u₁|, |d^{(2)}|, Rest |p_N − x₁ − d^{(2)}| (≈ d^{(3)}), Quotienten c₂ = |d^{(2)}|/ε², c₃ = |Rest|/ε³,
Test der Bilanz „eine Ableitung (N^{α}, α = 0.883) je Ordnung“: c₂ ~ N^α, c₃ ~ N^{2α}?
Daten: x₁ aus gemeinsam/lin_PiN_fK_bis_1000.json, p_N aus gemeinsam/pN_N100_bis_1000.json (Strang A)."""
import json, os, sys, time
import mpmath as mp
from flint import arb, acb, arb_poly, acb_poly, arb_mat, acb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = int(os.environ.get("PREC", "1000"))
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
K4 = d4["a"] if isinstance(d4, dict) else d4
lin = json.load(open(G + "lin_PiN_fK_bis_1000.json"))
PN = json.load(open(G + "pN_N100_bis_1000.json"))["p"]
mp.mp.dps = 60
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
aL = abs(Lm)
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 240) for j in range(240)] + [Lm, mp.conj(Lm)]
Kmp = [mp.mpf(c) for c in K4]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


fKw = [horner(Kmp, z) for z in kreis]


def supP(c):
    cm = [mp.mpf(x.mid().str(50, radius=False)) if hasattr(x, "mid") else mp.mpf(x) for x in c]
    return max(abs(horner(cm, z)) for z in kreis)


def coef(p, k):
    c = p.coeffs()
    return c[k] if k < len(c) else arb(0)


def tr(p, n):
    c = p.coeffs()[:n]
    return arb_poly(c) if c else arb_poly([0])


def inv_ser(p, n):
    y = arb_poly([1 / coef(p, 0)])
    m = 1
    while m < n:
        m = min(2 * m, n)
        y = tr(y * (arb_poly([2]) - tr(p, m) * y), m)
    return y


for N in [int(s) for s in sys.argv[1:]] or [100, 200, 300]:
    t0 = time.time()
    M = N + 2
    f = arb_poly([arb(c) for c in K[:M + 1]])
    pw = [arb_poly([1])]
    for _ in range(M):
        pw.append((pw[-1] * f).truncate(M))
    ex = [arb(1)]
    for k in range(1, M + 1):
        ex.append(ex[-1] / k)
    fp = tr(f.derivative(), M); ifp = inv_ser(fp, M)
    g = tr(arb_poly(ex[:M]) * ifp, M); gp = g.derivative(); fpp_f = tr(gp * ifp, M)
    Jm = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Jm[k, j] = coef(pw[j], k) + (coef(g, k - j) if j <= k else arb(0))
    Ji = Jm.inv()
    sol = lambda v: [x for x in (Ji * arb_mat([[coef(v, k)] for k in range(N + 1)])).entries()]
    tN = [Ji[k, N] for k in range(N + 1)]
    x1 = [arb(s) for s in lin[str(N)]["b"]]
    pN = [arb(s) for s in PN[str(N)]]
    # u₁ = x₁ − f als Keim bei 0 (mod x^M), u₁∘f = x₁∘f − eˣ
    u = tr(arb_poly(x1) - f, M)
    x1f = arb_poly([0])
    for j in range(N + 1):
        x1f += x1[j] * pw[j]
    uf = tr(x1f - arb_poly(ex[:M]), M)
    rho = coef(tr(uf + g * u, M), N + 1)
    # direkt
    duf = tr(tr(uf.derivative(), N + 1) * ifp, N + 1)                  # u'∘f
    uN = tr(u, N + 1)
    D2 = tr(2 * duf * uN + tr(fpp_f * tr(uN * uN, N + 1), N + 1), N + 1)
    d2 = [-x / 2 for x in sol(D2)]
    # zerlegt: v = u u'/(2f'), v∘f = (u∘f)(u'∘f)/(2g)
    v = tr(tr(u * u.derivative(), N + 1) * ifp, N + 1)
    v = arb_poly([c / 2 for c in v.coeffs()])
    vf = tr(tr(uf * duf, N + 1) * inv_ser(tr(g, N + 1), N + 1), N + 1)
    vf = arb_poly([c / 2 for c in vf.coeffs()])
    PiV = sol(tr(vf + g * v, N + 1))
    d2z = [PiV[k] - (N + 1) * rho * coef(u, 0) / (2 * coef(fp, 0)) * tN[k] for k in range(N + 1)]
    eps = max(abs(horner([mp.mpf(x.mid().str(50, radius=False)) for x in x1], z) - w) for z, w in zip(kreis, fKw))
    s2 = supP(d2); abw = supP([d2[k] - d2z[k] for k in range(N + 1)])
    kopf = supP([(N + 1) * rho * coef(u, 0) / (2 * coef(fp, 0)) * tN[k] for k in range(N + 1)])
    echt = [pN[k] - x1[k] for k in range(N + 1)]
    se = supP(echt)
    rest = supP([echt[k] - d2[k] for k in range(N + 1)])
    print(f"N={N} ({time.time()-t0:.0f}s): ε = {mp.nstr(eps,4)}, sup|p_N − x₁| = {mp.nstr(se,4)}, sup|d²| = {mp.nstr(s2,4)}"
          f" (Kopfanteil {mp.nstr(kopf,2)}, direkt↔zerlegt {mp.nstr(abw,2)}); Rest |p − x₁ − d²| = {mp.nstr(rest,4)};"
          f" c₂ = {mp.nstr(s2/eps**2,4)}, c₃ = {mp.nstr(rest/eps**3,4)}", flush=True)
