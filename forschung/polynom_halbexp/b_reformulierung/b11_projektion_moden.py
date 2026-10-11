"""b11: Lineares Modell der Auswahl. J-Projektion an f = f_K:  Π_N u = J_N^{-1} T_N[u∘f + (f'∘f)·u],
J_N δ = T_N[δ∘f + (f'∘f)·δ] auf Polynomen vom Grad ≤ N.  Π_N P = P für P ∈ P_N, Π_N k = 0 für k ∈ ker J_∞.
Gemessen: Modus-Inhalt von v = Π_N u − u in Kneser-Koordinaten, θ_v(w) = v(T(w)) · A'(f_K(T(w))),
Fourier-Koeffizienten bei e^{±2πiw} auf der Strecke w ∈ [−1.55, −0.55] (y = 0).
Testfunktionen u_s(z) = (1 − z/L)^s + (1 − z/L̄)^s (reell), s = Exponent der Singularität bei L, L̄.
Hypothese: Koeffizient ~ N^{a(s)};  bei p_N selbst (nichtlinear) wurde c_N ~ N^{1−2γ̄} gemessen (b10)."""
import json, sys, time
import mpmath as mp
from flint import arb, arb_poly, arb_mat, acb, ctx
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 700
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
Ac = [mp.mpf(c) for c in json.load(open(G + "kneser_slog_koeff_300.json"))["A"]]
dAc = [k * Ac[k] for k in range(1, len(Ac))]
Kmp = [mp.mpf(c) for c in K]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
gam = 2j * mp.pi / L

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

def T(w, z0):
    z = mp.mpc(z0)
    for _ in range(80):
        dz = (horner(Ac, z) - w) / horner(dAc, z)
        z -= dz
        if abs(dz) < mp.mpf(10) ** -50:
            return z
    raise RuntimeError

M, w0 = 64, mp.mpf("-1.55")
pts, z = [], mp.mpf(-0.5)
for j in range(M):
    w = w0 + mp.mpf(j) / M
    z = T(w, z)
    pts.append((w, mp.re(z)))
gewicht = [horner(dAc, horner(Kmp, z)) for w, z in pts]          # A'(f_K(z)) = 1/T'(w+1/2)
ex = {n: [mp.exp(-2j * mp.pi * n * w) for w, _ in pts] for n in (-1, 1)}


def baue(N, J=900):
    f = arb_poly([arb(c) for c in K[:N + 1]])
    pw = [arb_poly([1])]
    for _ in range(J):
        pw.append((pw[-1] * f).truncate(N + 1))
    ex_ser = [arb(1)]
    for k in range(1, N + 1):
        ex_ser.append(ex_ser[-1] / k)
    fp = arb_poly([arb(K[k]) * k for k in range(1, N + 2)])
    g = (arb_poly(ex_ser) * fp.inv_series(N + 1)).truncate(N + 1) if hasattr(fp, "inv_series") else None
    if g is None:   # Division von Potenzreihen von Hand: g = e^x / f'(x)
        c = [fp[k] for k in range(N + 1)]
        gk = []
        for n in range(N + 1):
            s = ex_ser[n] - sum((c[k] * gk[n - k] for k in range(1, n + 1)), arb(0))
            gk.append(s / c[0])
        g = arb_poly(gk)
    Jm = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Jm[k, j] = pw[j][k] + (g[k - j] if j <= k else 0)
    return pw, g, Jm.inv()


def proj(N, pw, g, Ji, uj):
    """uj: Taylor-Koeffizienten von u (Liste mpmath, Länge >= len(pw)); gibt Koeffizienten von Π_N u (mpmath)."""
    J = len(pw) - 1
    rhs = [arb(0)] * (N + 1)
    for j in range(J + 1):
        cj = arb(mp.nstr(uj[j], 70))
        for k in range(N + 1):
            rhs[k] += cj * pw[j][k]
    gu = (g * arb_poly([arb(mp.nstr(uj[j], 70)) for j in range(N + 1)])).truncate(N + 1)
    b = arb_mat([[rhs[k] + gu[k]] for k in range(N + 1)])
    x = Ji * b
    return [mp.mpf(x[k, 0].mid().str(70, radius=False)) for k in range(N + 1)]


def moden(N, coef, ufun):
    th = [(horner(coef, z) - ufun(z)) * gw for (w, z), gw in zip(pts, gewicht)]
    return {n: mp.fsum(t * e for t, e in zip(th, ex[n])) / M for n in (-1, 1)}


tests = {}
for s in (mp.mpf("1.5"), mp.mpf("2.5"), mp.mpf("3.5"), 1 + gam, mp.mpf("5.5")):
    def uj_fun(s=s, J=900):
        out, cL, cLb = [], mp.mpc(1), mp.mpc(1)
        for j in range(J + 1):
            out.append(mp.re(cL + cLb))
            cL = cL * (s - j) / (j + 1) * (-1 / L)
            cLb = cLb * (mp.conj(s) - j) / (j + 1) * (-1 / mp.conj(L))
        return out
    tests[mp.nstr(s, 5)] = (uj_fun(), lambda z, s=s: mp.re((1 - z / L) ** s + (1 - z / mp.conj(L)) ** mp.conj(s)))

res = {}
for N in (30, 50, 70, 90, 110):
    t0 = time.time()
    pw, g, Ji = baue(N)
    zeile = []
    for name, (uj, uf) in tests.items():
        co = proj(N, pw, g, Ji, uj)
        m = moden(N, co, uf)
        res.setdefault(name, {})[N] = m
        zeile.append(f"s={name}: c(-1)={mp.nstr(m[-1], 4)}")
    print(f"N={N} ({time.time()-t0:.0f}s): " + " | ".join(zeile), flush=True)

print("\nlokale Exponenten von |c(-1)| (zwischen aufeinanderfolgenden N):")
Ns = sorted(next(iter(res.values())).keys())
for name, d in res.items():
    ex_ = [mp.nstr(mp.log(abs(d[b][-1]) / abs(d[a][-1])) / mp.log(mp.mpf(b) / a), 4) for a, b in zip(Ns, Ns[1:])]
    print(f"  s={name}: {ex_}")
