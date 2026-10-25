# Fatou-Koordinaten und Hornabbildungen von g(u) = e^u - 1 (konjugiert zu eta^z bei z = e: z = e + e*u).
# phi(u) ~ -2/u + (1/3) log(u) + sum_{k>=1} c_k u^k  (formale Fatou-Koordinate, Konstante 0)
# phi_att: log(-u) (reell auf u<0), phi_rep: log(u) (reell auf u>0).
# Hornabbildung h(W) = phi_att(Psi_rep(W)); oben: h = W - i*pi/3 + sum a_n e^{2 pi i n W}.
from fractions import Fraction as F
import mpmath as mp

def formal_coeffs(K):
    """c_1..c_K der formalen Fatou-Koordinate (exakt, Brueche)."""
    N = K + 3  # Reihenlaenge
    # g(u)/u = (e^u-1)/u = sum u^j/(j+1)!
    fac = [F(1)]
    for j in range(1, N + 3):
        fac.append(fac[-1] * j)
    q = [F(1, 1) / fac[j + 1] for j in range(N + 1)]  # g/u
    def mul(a, b):
        r = [F(0)] * (N + 1)
        for i, x in enumerate(a):
            if x == 0: continue
            for j in range(N + 1 - i):
                r[i + j] += x * b[j]
        return r
    def inv(a):
        r = [F(0)] * (N + 1); r[0] = 1 / a[0]
        for n in range(1, N + 1):
            s = sum(a[k] * r[n - k] for k in range(1, n + 1))
            r[n] = -s / a[0]
        return r
    def logser(a):  # log(a), a[0]=1
        # log a = int a'/a
        da = [a[k + 1] * (k + 1) for k in range(N)] + [F(0)]
        t = mul(da, inv(a))
        r = [F(0)] * (N + 1)
        for k in range(N):
            r[k + 1] = t[k] / (k + 1)
        return r
    iq = inv(q)            # u/g
    lq = logser(q)         # log(g/u)
    # Glieder als Reihe in u (Index = Potenz):
    # -2/g + 2/u = (2/u)(1 - u/g) = 2 * (1 - iq)/u
    base = [F(0)] * (N + 1)
    t = [(-iq[k] if k else 1 - iq[0]) for k in range(N + 1)]
    for k in range(N):
        base[k] += 2 * t[k + 1]
    for k in range(N + 1):
        base[k] += F(1, 3) * lq[k]
    base[0] -= 1
    # g^k als Reihen
    g = [F(0)] + q[:N]  # g = u*q
    c = [F(0)] * (K + 1)
    gp = [g]  # gp[k-1] = g^k
    for k in range(2, K + 1):
        gp.append(mul(gp[-1], g))
    total = base[:]
    for k in range(1, K + 1):
        # Koeffizient von u^{k+1} in total + c_k (g^k - u^k) muss 0 sein; g^k - u^k = (k/2) u^{k+1} + ...
        assert total[k] == 0 or k == 0, (k, total[k])
        ck = -total[k + 1] / F(k, 2)
        c[k] = ck
        for j in range(N + 1):
            total[j] += ck * gp[k - 1][j]
        total[k] -= ck
    return c[1:]

_C = {}
def coeffs_mp(K):
    if K not in _C:
        _C[K] = [mp.mpf(x.numerator) / x.denominator for x in formal_coeffs(K)]
    return _C[K]

def phi_formal(u, K, att):
    c = coeffs_mp(K)
    s = mp.mpf(0); p = u
    for ck in c:
        s += ck * p; p *= u
    lg = mp.log(-u) if att else mp.log(u)
    return -2 / u + lg / 3 + s

def dphi_formal(u, K):
    c = coeffs_mp(K)
    s = mp.mpf(0); p = mp.mpf(1)
    for k, ck in enumerate(c, 1):
        s += k * ck * p; p *= u
    return 2 / u**2 + 1 / (3 * u) + s

def g(u):
    return mp.expm1(u)

def phi_att(u, K=30, delta=0.05, nmax=200000):
    n = 0
    while not (abs(u) < delta and mp.re(u) < 0 and abs(mp.im(u)) < -mp.re(u)):
        u = g(u); n += 1
        if n > nmax:
            raise RuntimeError('phi_att: keine Konvergenz')
    return phi_formal(u, K, True) - n

def psi_rep(W, K=30, R=40):
    """Umkehrung der abstossenden Fatou-Koordinate, ganze Funktion."""
    n = 0
    while mp.re(W - n) > -R:
        n += 1
    Wp = W - n
    u = -2 / Wp
    for _ in range(100):
        du = (phi_formal(u, K, False) - Wp) / dphi_formal(u, K)
        u -= du
        if abs(du) < mp.mpf(10) ** (-mp.mp.dps + 3) * abs(u):
            break
    for _ in range(n):
        u = g(u)
    return u

def horn(W, **kw):
    return phi_att(psi_rep(W, **{k: v for k, v in kw.items() if k in ('K', 'R')}),
                   **{k: v for k, v in kw.items() if k in ('K', 'delta')})
