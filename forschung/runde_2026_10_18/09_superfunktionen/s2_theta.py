"""theta = h(S(z)) - z auf waagrechten Linien, Fourier-Koeffizienten a_n, Vorhersage aus der Puiseux-Form.

Aufruf: python s2_theta.py <abb> [M] [y1] [y2]
"""
import numpy as np, sys, json, math, cmath
from scipy.special import gamma
from zoo import *
from s1_loesen import ABB

nm = sys.argv[1]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 8192
ys = [float(sys.argv[3]) if len(sys.argv) > 3 else 0.004, float(sys.argv[4]) if len(sys.argv) > 4 else 0.008]
a = ABB[nm][0]()
sol = {e: Loesung.laden(a, f'{nm}_{e}') for e in 'AB'}


def theta_linie(s, y, M):
    x = -0.5 + np.arange(M + 1) / M
    z = x + 1j * y
    S = s.S(z)
    ls, n = a.logsigma(S)
    im = np.unwrap(ls.imag)
    th = (ls.real + 1j * im) / a.loglam - z
    return x, th


def ableitung(s, z0, r=0.05, K=16):
    w = np.exp(2j * np.pi * np.arange(K) / K)
    v = s.Sweit(z0 + r * w)
    return (v * np.conj(w)).mean() / r


out = {'abb': nm, 'L': [a.L.real, a.L.imag], 'lam': [a.lam.real, a.lam.imag], 'M': M, 'y': ys}
koeff = {}
for e in 'AB':
    s = sol[e]
    for y in ys:
        x, th = theta_linie(s, y, M)
        per = th[-1] - th[0]
        c = np.fft.fft(th[:-1]) / M * np.exp(-2j * np.pi * np.arange(M) * (-0.5))  # Bezug x=0
        n = np.arange(M); n[n > M // 2] -= M
        an = c * np.exp(2 * np.pi * n * y)
        koeff[(e, y)] = (n, an, c)
        neg = np.max(abs(c[(n < 0) & (n > -M // 8)]))
        print(f'{e} y={y}: Periodenfehler |th(1/2)-th(-1/2)| = {abs(per):.1e},  max|roh c_-m| = {neg:.1e}', flush=True)

n, aA, _ = koeff[('A', ys[0])]
_, aB, _ = koeff[('B', ys[0])]
_, aA2, _ = koeff[('A', ys[1])]
pos = lambda arr: arr[(n >= 0)][:M // 2]
aA, aB, aA2 = pos(aA), pos(aB), pos(aA2)
nn = np.arange(len(aA))
print('a_0 =', aA[0], ' a_1 =', aA[1], ' a_2 =', aA[2])

# Vorhersage
s = sol['A']
kr = a.krit
z1 = -1.0
A1 = ableitung(s, z1)  # S'(z1), S(z1) = v
print('S(z1) =', s.Sweit(np.array([z1 + 0j]))[0], ' A = S\'(z1) =', A1)
res = {}
if math.isinf(kr['d']):
    c1 = cmath.log(A1.real / (2 * np.pi)) + 1j * np.pi / 2
    lnb = a.lnb
    def vorh(k):
        return -a.hp(np.array([(c1 - math.log(k)) / lnb]))[0] / (k * lnb)
    print('c =', c1)
else:
    d = kr['d']; hpc = a.hp(np.array([kr['wc'] + 0j]))[0]
    C1 = hpc * (A1.real * 1j / (2 * np.pi * kr['a'])) ** (1.0 / d) / gamma(-1.0 / d)
    print("h'(w_c) =", hpc, '  C1 =', C1)
    def vorh(k):
        return C1 * k ** (-1 - 1.0 / d)
    res['C1'] = [C1.real, C1.imag]; res['hp_wc'] = [hpc.real, hpc.imag]
res['A'] = A1.real
tab = []
ks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
ks = [k for k in ks if k < len(aA)]
print(' n      |a_n|        lokExp   arg a_n     a_n/vorh                   |A-B|/|a|   |y1-y2|/|a|')
for k in ks:
    le = -math.log(abs(aA[2 * k]) / abs(aA[k])) / math.log(2) if 2 * k < len(aA) else float('nan')
    q = aA[k] / vorh(k)
    dAB = abs(aA[k] - aB[k]) / abs(aA[k]); dy = abs(aA[k] - aA2[k]) / abs(aA[k])
    print(f'{k:5d}  {abs(aA[k]):.6e}  {le:.4f}  {np.angle(aA[k]):+.5f}  {q.real:.7f}{q.imag:+.7f}i   {dAB:.1e}   {dy:.1e}')
    tab.append(dict(n=k, a=[aA[k].real, aA[k].imag], quot=[q.real, q.imag], dAB=dAB, dy=dy, lokexp=le))
res['tab'] = tab
res['a'] = [[v.real, v.imag] for v in aA[:1025]]
out.update(res)
json.dump(out, open(os.path.join(DATEN, f'theta_{nm}.json'), 'w'), indent=1)
