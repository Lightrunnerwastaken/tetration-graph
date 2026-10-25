# Kanonische (Fatou-/Ecalle-) Fortsetzung von ausgewogener Tetration und linker Pentation.
# Beide sind in Log-Koordinaten die Iteration von F(m) = m e^m (Umkehrung: Lambert W).
#   ausgewogen:  a[4]b   = exp( F^[log2 b](ln a) )
#   links:       x[5']t  = exp( F^[t](x ln x) / x )
#   ausgewogene Pentation-Schritt: P(y) = y[4]y,  ln P(e^l) = F^[l/ln2](l)
# F^[t] ueber die Fatou-Koordinate Phi(m) = -1/m + (1/2) ln|m| + sum c_k m^k (fatou_reihe.py),
# m > 0: abstossende Seite (Rueckwaerts mit W), m < 0: anziehende Seite (vorwaerts mit F).
# Aufruf: python kanonisch.py <dps> <n> <K>   (zwei Einstellungen = Gegenprobe)
import sys
from fractions import Fraction as Q
import mpmath as mp

dps = int(sys.argv[1]); n = int(sys.argv[2]); K = int(sys.argv[3])
mp.mp.dps = dps

C = [Q(-1, 6), Q(1, 16), Q(-19, 540), Q(1, 48), Q(-41, 4200), Q(37, 103680), Q(18349, 3175200),
     Q(-443, 80640), Q(-55721, 21555072), Q(84317, 6912000), Q(-2594833561, 399567168000),
     Q(-152043613, 5748019200), Q(830066563, 17346700800), Q(209959420499, 3905587445760),
     Q(-8125100271473651, 30011126745600000), Q(-3290126451127, 502146957312000),
     Q(35263493868569, 21662202691584), Q(-82344139152648353, 51218989645824000)]
C = [mp.mpf(c.numerator)/c.denominator for c in C[:K]]

def phi_as(m):
    r = -1/m + mp.log(abs(m))/2
    p = m
    for c in C:
        r += c*p; p *= m
    return r

def dphi_as(m):
    r = 1/m**2 + 1/(2*m)
    p = mp.mpf(1)
    for k, c in enumerate(C, 1):
        r += k*c*p; p *= m
    return r

def inv_as(v, m0):
    m = m0
    for _ in range(60):
        dm = (phi_as(m) - v)/dphi_as(m)
        m -= dm
        if abs(dm) < mp.mpf(10)**(-dps + 3)*abs(m):
            break
    return m

F = lambda m: m*mp.exp(m)
W = lambda m: mp.lambertw(m).real

def Fit(t, m):
    m = mp.mpf(m)
    if m == 0:
        return mp.mpf(0)
    if m > 0:
        mu = m
        for _ in range(n): mu = W(mu)
        mu2 = inv_as(phi_as(mu) + t, mu)
        for _ in range(n): mu2 = F(mu2)
        return mu2
    else:
        mu = m
        for _ in range(n): mu = F(mu)
        mu2 = inv_as(phi_as(mu) + t, mu)
        for _ in range(n): mu2 = W(mu2)
        return mu2

def bal4(a, b):
    return mp.exp(Fit(mp.log(b)/mp.log(2), mp.log(a)))

def L5(x, t):
    x = mp.mpf(x)
    return mp.exp(Fit(t, x*mp.log(x))/x)

out = []
def pr(name, v):
    out.append((name, v)); print(name, mp.nstr(v, dps - 8))

pr('Probe bal4(2,2)=4', bal4(2, 2))
pr('Probe bal4(2,4)=256', bal4(2, 4))
pr('Probe L5(2,1)=16', L5(2, 1))
pr('Probe L5(1/2,1)=(1/2)^((1/2)^(1/2))', L5(mp.mpf(1)/2, 1))
pr('bal4(2,3)', bal4(2, 3))
pr('bal4(2,sqrt2)', bal4(2, mp.sqrt(2)))
pr('bal4(e,3)', bal4(mp.e, 3))
pr('bal4(1/2,3)', bal4(mp.mpf(1)/2, 3))
pr('L5(2,1/2)', L5(2, mp.mpf(1)/2))
pr('L5(2,3/2)', L5(2, mp.mpf(3)/2))
pr('L5(e,1/2)', L5(mp.e, mp.mpf(1)/2))
pr('L5(1/2,1/2)', L5(mp.mpf(1)/2, mp.mpf(1)/2))
pr('F^[1/2](1)', Fit(mp.mpf(1)/2, 1))
pr('F^[1/2](-1/2)', Fit(mp.mpf(1)/2, -mp.mpf(1)/2))
# Halbgruppe: F^[1/2](F^[1/2](1)) = F(1) = e
pr('Gruppenprobe F^[1/2]^2(1) - e', Fit(mp.mpf(1)/2, Fit(mp.mpf(1)/2, 1)) - mp.e)
# ausgewogene Pentation: (ln P(e^l) - l)/l^3 -> 1/ln 2
for l in ('0.02', '0.005', '-0.005', '-0.02'):
    l = mp.mpf(l)
    v = Fit(l/mp.log(2), l)
    pr('ausgew. Pentation (lnP-l)/l^3 bei l=%s' % mp.nstr(l, 3), (v - l)/l**3)
pr('1/ln2', 1/mp.log(2))
