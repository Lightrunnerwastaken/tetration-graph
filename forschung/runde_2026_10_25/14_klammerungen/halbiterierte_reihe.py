# Formale halbe Iterierte g = F^[1/2] von F(m) = m e^m: g(g(m)) = F(m), g = m + m^2/2 + ...
# (fuer tangential-identische Keime eindeutig). Waechst |g_k| faktoriell, so gibt es KEINE an 0
# analytische Halbiterierte; dann ist die ausgewogene Tetration a[4]b fuer b != 2^n bei a = 1
# zwar C^unendlich, aber nicht analytisch (beidseitige Fatou-Fortsetzungen kleben nicht analytisch).
from fractions import Fraction as Q
import math, sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30

def mul(a, b):
    c = [Q(0)]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]:
                    c[i+j] += x*b[j]
    return c

def compose(f, g):
    # f(g), g[0] = 0
    r = [Q(0)]*N
    p = [Q(1)] + [Q(0)]*(N - 1)
    for k in range(1, N):
        p = mul(p, g)
        if f[k]:
            r = [x + f[k]*y for x, y in zip(r, p)]
    return r

Fc = [Q(0)] + [Q(1, math.factorial(k - 1)) for k in range(1, N)]   # m e^m = sum m^k/(k-1)!
g = [Q(0)]*N; g[1] = Q(1)
for k in range(2, N):
    h = compose(g, g)
    g[k] = (Fc[k] - h[k])/2
for k in range(2, N):
    v = float(g[k])
    norm = abs(v)*(2*math.pi)**k/math.factorial(k - 2) if v else 0.0
    print(k, g[k] if k < 12 else '', '%.6e' % v, ' |g_k|(2pi)^k/(k-2)! = %.4f' % norm)
