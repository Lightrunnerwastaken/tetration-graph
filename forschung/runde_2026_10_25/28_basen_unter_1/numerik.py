"""Zahlen fuer BERICHT.md (Gebiet 28).  Aufruf: python numerik.py teil   (teil = werte | sigma | kritisch | rand | realteil | alle)

werte    : T_cos(1/4), T_cos(3/4), T_cos'(0) fuer fuenf Basen, zwei Einstellungen (dps 40/tol 1e-3 gegen dps 60/tol 1e-5)
sigma    : sigma auf R streng wachsend, Bild (sigma(1)/lambda^3, sigma(1)/lambda^2)
kritisch : kritische Punkte von T_cos bei n - delta, delta = arctan(|ln|lambda||/pi)/pi
rand     : b -> e^{-e}+ : Verhalten von T_cos(1/4)
realteil : Re T_reg ist keine Tetration
"""
import sys
import mpmath as mp
import schroeder as S


def einstellung(dps, tol):
    mp.mp.dps = dps + 10
    S.NSER = 40 if dps <= 40 else 60
    return mp.mpf(tol)


def Tcos(P, x, tol):
    s1 = S.sigma(P, mp.mpf(1), tol)
    return S.psi(P, abs(P['lam']) ** x * s1 * mp.cospi(x), tol)


def teil_werte():
    print('== werte ==')
    for b in ['1/10', '3/10', '1/2', '7/10', '9/10']:
        res = []
        for dps, tol in [(40, '1e-3'), (60, '1e-5')]:
            t = einstellung(dps, tol)
            P = S.setup(mp.mpf(int(b.split('/')[0])) / int(b.split('/')[1]))
            v = [Tcos(P, mp.mpf(1) / 4, t), Tcos(P, mp.mpf(3) / 4, t),
                 mp.diff(lambda x: Tcos(P, x, t), 0)]
            res.append((P, v))
        (P, v), (P2, v2) = res
        mp.mp.dps = 50
        d = max(abs(x - y) for x, y in zip(v, v2))
        print(f'b={b}  a={mp.nstr(P2["a"], 32)}  lam={mp.nstr(P2["lam"], 32)}')
        print(f'   T(1/4)={mp.nstr(v2[0], 32)}  T(3/4)={mp.nstr(v2[1], 32)}  T\'(0)={mp.nstr(v2[2], 25)}   Abweichung zweier Einstellungen {mp.nstr(d, 3)}')


def teil_sigma():
    print('== sigma ==')
    t = einstellung(30, '1e-3')
    for b in [mp.mpf(1) / 10, mp.mpf(1) / 2, mp.mpf(9) / 10]:
        P = S.setup(b)
        s1 = S.sigma(P, mp.mpf(1), t)
        xs = [mp.mpf(k) / 8 for k in range(-15, 49)]
        vs = [S.sigma(P, x, t) for x in xs]
        mono = all(vs[i] < vs[i + 1] for i in range(len(vs) - 1))
        print(f'b={mp.nstr(b,3)}: wachsend auf [-15/8, 6] (Schritt 1/8): {mono};  sigma(6)={mp.nstr(vs[-1], 15)}  sigma(1)/lam^2={mp.nstr(s1 / P["lam"]**2, 15)};'
              f'  sigma(-15/8)={mp.nstr(vs[0], 15)}  sigma(1)/lam^3={mp.nstr(s1 / P["lam"]**3, 15)}')


def teil_kritisch():
    print('== kritisch ==')
    t = einstellung(30, '1e-3')
    for b in [mp.mpf(1) / 10, mp.mpf(1) / 2, mp.mpf(9) / 10]:
        P = S.setup(b)
        delta = mp.atan(abs(mp.log(abs(P['lam']))) / mp.pi) / mp.pi
        g = lambda x: mp.diff(lambda y: Tcos(P, y, t), x)
        r = [mp.findroot(g, n - delta + mp.mpf('0.01')) for n in (0, 1, 2)]
        print(f'b={mp.nstr(b,3)}: delta={mp.nstr(delta, 20)}  Nullstellen von T\' bei {[mp.nstr(x, 20) for x in r]}'
              f'  T(-delta)={mp.nstr(Tcos(P, -delta, t), 15)} > T(0)=1')


def teil_rand():
    print('== rand b -> e^{-e}+ ==')
    t = einstellung(25, '1e-4')
    S.NSER = 60
    for tt in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['1e-1', '3e-2', '1e-2', '3e-3', '1e-3']):
        b = mp.exp(-mp.e) * (1 + mp.mpf(tt))
        P = S.setup(b)
        s1 = S.sigma(P, mp.mpf(1), t)
        vals = [Tcos(P, mp.mpf(x), t) for x in ('0.25', '0.75', '1.25', '2.25')]
        print(f't={tt}: lam={mp.nstr(P["lam"], 12)} sigma(1)={mp.nstr(s1, 12)}  T(1/4)={mp.nstr(vals[0], 14)} T(3/4)={mp.nstr(vals[1], 14)}'
              f' T(5/4)={mp.nstr(vals[2], 14)} T(9/4)={mp.nstr(vals[3], 14)}', flush=True)


def teil_realteil():
    print('== realteil ==')
    t = einstellung(30, '1e-3')
    P = S.setup(mp.mpf(1) / 2)
    for x in [mp.mpf(1) / 4, mp.mpf(1) / 2]:
        r0 = S.T_reg(P, x).real
        r1 = S.T_reg(P, x + 1).real
        print(f'x={mp.nstr(x,3)}: Re T_reg(x+1) - b^(Re T_reg(x)) = {mp.nstr(r1 - S.f(P, r0), 12)}')


if __name__ == '__main__':
    teil = sys.argv[1] if len(sys.argv) > 1 else 'alle'
    for name, fn in [('werte', teil_werte), ('sigma', teil_sigma), ('kritisch', teil_kritisch),
                     ('realteil', teil_realteil), ('rand', teil_rand)]:
        if teil in (name, 'alle'):
            fn()
