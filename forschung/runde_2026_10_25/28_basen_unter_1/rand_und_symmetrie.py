"""(1) Symmetrie T_cos(-x) = g^{-x}(T_cos(x)), g = f o f mit reellem regulaerem Fluss g^s = Psi(lambda^{2s} sigma(.)).
(2) Randverhalten b = e^{-e}(1+t):  sigma_b(1) / sqrt(t)  und  (T_cos(x) - a_b) / (sqrt(t) cot(pi x))  gegen sqrt(6) e^{-3/2}.
Aufruf: python rand_und_symmetrie.py [t-Liste, Komma]
"""
import sys
import mpmath as mp
import schroeder as S

mp.mp.dps = 30
S.NSER = 60
tol = mp.mpf('1e-4')


def Tcos(P, x):
    s1 = S.sigma(P, mp.mpf(1), tol)
    return S.psi(P, abs(P['lam']) ** x * s1 * mp.cospi(x), tol)


def gflow(P, s, z):
    return S.psi(P, abs(P['lam']) ** (2 * s) * S.sigma(P, z, tol), tol)


print('== Symmetrie ==')
for b in [mp.mpf(1) / 10, mp.mpf(1) / 2, mp.mpf(9) / 10]:
    P = S.setup(b)
    for x in [mp.mpf('0.3'), mp.mpf('0.77'), mp.mpf('1.4')]:
        lhs = Tcos(P, -x)
        rhs = gflow(P, -x, Tcos(P, x))
        print(f'b={mp.nstr(b,3)} x={mp.nstr(x,3)}: T(-x)={mp.nstr(lhs, 20)}  g^(-x)(T(x))={mp.nstr(rhs, 20)}  diff={mp.nstr(abs(lhs-rhs), 3)}')

print('== Rand ==  kappa = sqrt(6) e^{-3/2} =', mp.nstr(mp.sqrt(6) * mp.exp(-1.5), 15))
ts = sys.argv[1].split(',') if len(sys.argv) > 1 else ['1e-2', '1e-3', '3e-4']
for tt in ts:
    t = mp.mpf(tt)
    P = S.setup(mp.exp(-mp.e) * (1 + t))
    s1 = S.sigma(P, mp.mpf(1), tol)
    row = [f'sigma(1)/sqrt(t)={mp.nstr(s1 / mp.sqrt(t), 10)}']
    for x in ['0.25', '0.4', '0.7']:
        x = mp.mpf(x)
        q = (Tcos(P, x) - P['a']) / (mp.sqrt(t) * mp.cot(mp.pi * x))
        row.append(f'x={mp.nstr(x,2)}: {mp.nstr(q, 10)}')
    print(f't={tt}: ' + '  '.join(row), flush=True)
