"""Die Identitäten für eine allgemeine Basis b > η (hier b = 2; b = e als Gegenprobe).

λ = L ln b,  σ(b^w) = λ σ(w),  h(w) = n + log σ(u_n)/log λ,  u_{k+1} = log_b u_k,  λ^{h(w)} = σ(w),
θ(x) = h(sexp x) − x.  Geprüft:
  Z  M[θ²] = M[θ]²,  M[λ^θ]·M[λ^{−θ}] = 1,  M[log(1+θ')] = 0,  M[1/(1+θ')] = 1,  M[θ e^{2πix}] = 0
  V1 ∫_ℝ [b^{sexp(it)} − sexp(it)] dt = 2 Im L
  V3 ∫_0^1 [ln sexp' − log λ · h(sexp)] dx ≡ log log λ − i ln b ∫_0^∞ (sexp(it) − L) dt   (mod 2πi)
Aufruf: python s5_basis_b.py BASIS DPS     (BASIS: 2 oder e; nutzt Motor für Taylor-Daten, θ aus etappe2)
"""
import json
import sys
import time

import mpmath as mp

from k3 import HIER, kern

arg_b = sys.argv[1] if len(sys.argv) > 1 else "2"
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
mp.mp.dps = dps + 10
basis = "exp(1)" if arg_b == "e" else arg_b
b = mp.e if arg_b == "e" else mp.mpf(arg_b)
lb = mp.log(b)
L = -mp.lambertw(-lb) / lb
if mp.im(L) < 0:
    L = mp.conj(L)
lam = L * lb
llam = mp.log(lam)
K = 70

# ---- Taylor-Daten (Kreise r = 0.2 um −0.4 … 0.4)
tf = HIER / "daten" / f"taylor_b{arg_b}_{dps + 30}_256.json"
if not tf.is_file():
    M = 256; r = mp.mpf("0.2"); out = {}
    for k in range(-2, 3):
        c = mp.mpf(k) / 5
        zs = [c + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
        v = kern.werte(basis, [f"sexp({kern.gpzahl(z, dps + 40)})" for z in zs], dps + 30)
        out[mp.nstr(c, 5)] = [mp.nstr(mp.re(mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * n) / M) for j in range(M)) / M / r ** n), dps + 35) for n in range(M // 2)]
    tf.write_text(json.dumps(out), encoding="utf-8")
ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in json.loads(tf.read_text(encoding="utf-8")).items()}

# ---- Schröder σ an L und Ψ
g = [mp.mpc(0)] + [L * lb ** k / mp.factorial(k) for k in range(1, K + 1)]


def mul(a, c):
    o = [mp.mpc(0)] * (K + 1)
    for i in range(K + 1):
        if a[i] != 0:
            for j in range(K + 1 - i):
                o[i + j] += a[i] * c[j]
    return o


gp = [None, g]
for _ in range(2, K + 1):
    gp.append(mul(gp[-1], g))
S = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
for m in range(2, K + 1):
    S[m] = mp.fsum(S[k] * gp[k][m] for k in range(1, m)) / (lam - lam ** m)
dS = [k * S[k] for k in range(1, K + 1)]
P = [L, mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
for k in range(2, K + 1):
    q = [mp.mpc(0)] + [lb * x for x in P[1:k]] + [mp.mpc(0)] * (K + 1 - k)
    e = [mp.mpc(1)] + [mp.mpc(0)] * k
    for m in range(1, k + 1):
        e[m] = mp.fsum(j * q[j] * e[m - j] for j in range(1, m + 1)) / m
    P[k] = L * e[k] / (lam ** k - lam)


def Psi(u):
    n = 0
    while abs(u) > mp.mpf("0.02"):
        u /= lam; n += 1
    w = mp.polyval(P[::-1], u)
    for _ in range(n):
        w = b ** w
    return w


def h_d(w):
    """h(w) und h'(w)."""
    n = 0; prod = mp.mpc(1)
    tol = mp.mpf(10) ** (-mp.mp.dps)
    while abs(w - L) > mp.mpf("1e-3") or abs(w - L) ** (K + 1) > tol * abs(w - L):
        prod /= (w * lb)
        w = mp.log(w) / lb; n += 1
    u = w - L
    su = mp.polyval(S[::-1], u)
    return n + mp.log(su) / llam, prod * mp.polyval(dS[::-1], u) / su / llam


def a_lesen():
    d = json.loads((HIER.parent / "etappe2" / "daten" / f"a_{arg_b}_80_0.003_16384.json").read_text(encoding="utf-8"))
    return [mp.mpc(mp.mpf(x[0]), mp.mpf(x[1])) for x in d["a"]]


A = a_lesen()


def sexp_oben(z):
    q = mp.expjpi(2 * z); th = mp.mpc(0); qn = mp.mpc(1)
    for ak in A:
        th += ak * qn; qn *= q
        if abs(qn) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
    return Psi(mp.exp(llam * (z + th)))


def sexp(z):
    if mp.im(z) <= mp.mpf("0.09"):
        c = min(ZEN, key=lambda c: abs(z - c))
        if abs(z - c) <= mp.mpf("0.1"):
            return mp.polyval(ZEN[c][::-1], z - c)
    return sexp_oben(z)


def reihe(x):
    c = min(ZEN, key=lambda c: abs(x - c)); co = ZEN[c]; u = x - c
    s1 = mp.polyval(co[:0:-1], u) * u if c == 0 else mp.polyval(co[::-1], u) - 1
    return s1, mp.polyval([k * co[k] for k in range(len(co) - 1, 0, -1)], u)


def entwirren(w, per):
    out = [w[0]]
    for v in w[1:]:
        mr = int(mp.nint(mp.re((out[-1] - v) / per))) if mp.re(per) != 0 else int(mp.nint(mp.im((out[-1] - v) / per)))
        best = min([mr + m for m in range(-3, 4)], key=lambda m: abs(v + m * per - out[-1]))
        out.append(v + best * per)
    return out


def nullmoden(k, tmax=5):
    hh = mp.mpf(2) ** (-k); n = int(tmax / hh)
    th, lg, W, X = [], [], [], []
    for j in range(-n, n + 1):
        t = j * hh; a = mp.pi * mp.sinh(t)
        dl = 1 / (1 + mp.exp(-a)); dr = 1 / (1 + mp.exp(a))
        W.append(hh * mp.pi * mp.cosh(t) * dl * dr); X.append(dl)
        xs = dl if dl <= mp.mpf("0.5") else -dr
        s1, ds = reihe(xs)
        ell = mp.log1p(s1) / lb
        hv, hd = h_d(ell)
        th.append(1 + hv - xs)
        lg.append(mp.log(ds * hd / ((1 + s1) * lb)))
    mid = len(th) // 2
    per = 2j * mp.pi / llam
    th = entwirren(th[mid::-1], per)[::-1][:-1] + entwirren(th[mid:], per)
    lg = entwirren(lg[mid::-1], 2j * mp.pi)[::-1][:-1] + entwirren(lg[mid:], 2j * mp.pi)
    Mf = lambda f: mp.fsum(w * f(i) for i, w in enumerate(W))
    mt = Mf(lambda i: th[i])
    return {"M_theta": mt, "a0": A[0], "var": Mf(lambda i: th[i] ** 2) - mt ** 2,
            "lam": Mf(lambda i: lam ** th[i]) * Mf(lambda i: lam ** (-th[i])) - 1,
            "log": Mf(lambda i: lg[i]), "inv": Mf(lambda i: mp.exp(-lg[i])) - 1,
            "mode+1": Mf(lambda i: th[i] * mp.expjpi(2 * X[i])),
            "kappa": Mf(lambda i: lam ** th[i])}


if __name__ == "__main__":
    t0 = time.time()
    print(f"Basis {arg_b}: L = {mp.nstr(L, 25)}, λ = {mp.nstr(lam, 20)}")
    for k in (6, 7):
        r = nullmoden(k)
        print(f"k={k}: M[θ] = {mp.nstr(r['M_theta'], 25)}  (a0 etappe2 {mp.nstr(r['a0'], 15)}; Diff/(2πi/logλ) = {mp.nstr((r['M_theta']-r['a0'])/(2j*mp.pi/llam), 6)})")
        for key in ("var", "lam", "log", "inv", "mode+1"):
            print(f"     {key:7s} |Abw| = {mp.nstr(abs(r[key]), 3)}")
        print(f"     κ_b = M[λ^θ] = {mp.nstr(r['kappa'], 30)}", flush=True)
    q = lambda f: mp.quad(f, [0, mp.mpf("0.09"), 0.5, 2, 6, 15, 35, 70, 140])
    v1 = q(lambda t: b ** sexp(mp.mpc(0, t)) - sexp(mp.mpc(0, t)))
    print(f"V1: ∫_ℝ[b^s − s]dt = {mp.nstr(2*mp.re(v1), 30)}  2 Im L = {mp.nstr(2*mp.im(L), 30)}  |Abw| = {mp.nstr(abs(2*mp.re(v1)-2*mp.im(L)), 3)}")
    K0 = q(lambda t: sexp(mp.mpc(0, t)) - L)
    # ∫_0^1 ln sexp' (glatt); [1/2,1] über ln sexp'(x) = ln ln b + ln b · sexp(x−1) + ln sexp'(x−1)
    def lnd(x):
        if x <= 0.5:
            return mp.log(reihe(mp.mpf(x))[1])
        s1, ds = reihe(mp.mpf(x) - 1)
        return mp.log(lb) + lb * (1 + s1) + mp.log(ds)
    lhs = mp.quad(lnd, [0, 0.25, 0.5, 0.75, 1]) - llam * (r["M_theta"] + mp.mpf(1) / 2)
    rhs = mp.log(llam) - 1j * lb * K0
    print(f"V3: (links − rechts)/(2πi) = {mp.nstr((lhs - rhs) / (2j * mp.pi), 25)}  [ganzzahlig erwartet]")
    print(f"Laufzeit {time.time()-t0:.0f}s")
