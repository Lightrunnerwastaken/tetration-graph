# Linke Hexation zur Basis 2:  2[6']0 = 2,  2[6'](n+1) = (2[6']n)[5'] 2.
# In Log-Koordinaten ist der Schritt  G(z) = ln(e^z [5'] 2) = F^[3](z) e^{-z}  (F(m) = m e^m, ganze Iterierte!),
# parabolisch: G(z) = z + 2 z^2 + 5 z^3 + ...  (a = x = 2, b = x^2 + x/2 = 5, beta = 1 - b/a^2 = -1/4).
# Kanonische Fortsetzung ueber die abstossende Fatou-Koordinate (Rueckwaerts mit G^{-1} per Newton).
# Aufruf: python hexation_links.py <dps> <n> <K>
import sys
from fractions import Fraction as Q
import math
import mpmath as mp

dps = int(sys.argv[1]); n = int(sys.argv[2]); K = int(sys.argv[3])
mp.mp.dps = dps
N = K + 3

def mul(a, b):
    c = [Q(0)]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]: c[i+j] += x*b[j]
    return c

def compose(f, g):
    r = [Q(0)]*N; p = [Q(1)] + [Q(0)]*(N - 1)
    for k in range(1, N):
        p = mul(p, g)
        if f[k]: r = [x + f[k]*y for x, y in zip(r, p)]
    return r

Fc = [Q(0)] + [Q(1, math.factorial(k - 1)) for k in range(1, N)]
emz = [Q((-1)**k, math.factorial(k)) for k in range(N)]
G = mul(compose(Fc, compose(Fc, Fc)), emz)
a, b = G[2], G[3]
beta = 1 - b/(a*a)
print('G Koeff.:', [str(x) for x in G[:6]], ' beta =', beta)

# Fatou-Reihe: Phi(z) = -1/(a z) + beta ln z + sum_{k=1}^K c_k z^k, Phi(G) - Phi = 1.
# u := G/z - 1 (Reihe ab z^1);  -1/(aG) + 1/(az) = (1/(a z)) (1 - 1/(1+u));  ln(G/z) = ln(1+u).
u = G[1:] + [Q(0)]          # G/z
u[0] -= 1
def inv1p(u):               # 1/(1+u)
    r = [Q(1)] + [Q(0)]*(N - 1); t = [Q(1)] + [Q(0)]*(N - 1)
    for k in range(1, N):
        t = mul(t, [-x for x in u]); r = [x + y for x, y in zip(r, t)]
    return r
def log1p(u):
    r = [Q(0)]*N; t = [Q(1)] + [Q(0)]*(N - 1)
    for k in range(1, N):
        t = mul(t, u); r = [x + Q((-1)**(k+1), k)*y for x, y in zip(r, t)]
    return r
one_minus = [-x for x in inv1p(u)]; one_minus[0] += 1     # 1 - 1/(1+u), ab z^1
base = [Q(0)]*N
for i in range(N - 1):
    base[i] = one_minus[i+1]/a      # (1/(a z)) * (...), Index = Potenz
lg = log1p(u)
base = [x + beta*y for x, y in zip(base, lg)]
base[0] -= 1
# sum_k c_k (G^k - z^k): Koeffizient bei z^{k+1} ist k a c_k + (aeltere)
powG = [[Q(1)] + [Q(0)]*(N - 1)]
for k in range(1, K + 1): powG.append(mul(powG[-1], G))
cs = []
for k in range(1, K + 1):
    tot = base[k + 1 - 0] if False else base[k]   # Ordnung z^k der Restgleichung
    # Beitrag schon bestimmter c_j zur Ordnung z^k
    for j, c in enumerate(cs, 1):
        d = list(powG[j]); d[j] -= 1
        tot += c*d[k]
    # neues c_{k-1}?  Ordnung z^k bestimmt c_{k-1} (Leitterm (k-1) a c_{k-1} z^k)
    if k == 1:
        assert tot == 0, tot   # Ordnung z^1 legt beta fest
        continue
    cnew = -tot/((k - 1)*a)
    cs.append(cnew)
print('c_k:', [str(c) for c in cs[:6]])
C = [mp.mpf(c.numerator)/c.denominator for c in cs]
A = mp.mpf(a.numerator)/a.denominator; B = mp.mpf(beta.numerator)/beta.denominator

def phi(z):
    r = -1/(A*z) + B*mp.log(z); p = z
    for c in C: r += c*p; p *= z
    return r
def Gm(z):
    w = z
    for _ in range(3): w = w*mp.exp(w)
    return w*mp.exp(-z)
def Ginv(v):
    # G ist auf [0, 1] streng wachsend, G(z) >= z, G(1) ~ 1e19: Wurzel in (0, min(v,1)]
    hi = min(v, mp.mpf(1))
    return mp.findroot(lambda z: Gm(z) - v, (mp.mpf(0), hi), solver='anderson', tol=mp.mpf(10)**(-2*dps))
def Git(t, z):
    mu = z
    for _ in range(n): mu = Ginv(mu)
    v = phi(mu) + t
    m2 = mu
    for _ in range(100):
        h = mp.mpf(10)**(-dps//2)
        d = (phi(m2 + h*m2) - phi(m2))/(h*m2)
        step = (phi(m2) - v)/d; m2 -= step
        if abs(step) < mp.mpf(10)**(-dps + 5)*m2: break
    for _ in range(n): m2 = Gm(m2)
    return m2

z0 = mp.log(2)
print('Probe ln(2[6\']1) = 1024 ln2:', mp.nstr(Git(1, z0), 20), mp.nstr(1024*mp.log(2), 20))
v = Git(mp.mpf(1)/2, z0)
print("2[6']1/2 =", mp.nstr(mp.exp(v), dps - 15), '  ln =', mp.nstr(v, dps - 15))
print('Halbgruppe:', mp.nstr(Git(mp.mpf(1)/2, v) - 1024*mp.log(2), 5))
