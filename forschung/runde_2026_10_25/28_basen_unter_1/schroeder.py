"""Reelle Tetrationen fuer e^{-e} < b < 1 ueber die Schroeder-Funktion am anziehenden Fixpunkt a.

f(z) = b^z, a = W0(c)/c mit c = -ln b, lambda = a ln b in (-1, 0).
sigma(f(z)) = lambda sigma(z), sigma'(a) = 1;  Psi = sigma^{-1}.
Jede stetige reelle Tetration:  T_p(x) = Psi(|lambda|^x p(x)),  p(x+1) = -p(x), p(0) = sigma(1).
Kosinus-Wahl:  p(x) = sigma(1) cos(pi x)  ->  T_cos.

Aufruf: python schroeder.py [dps]
"""
import sys
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 40
mp.mp.dps = DPS + 10
NSER = 40  # Laenge der Schroeder-Reihe


def setup(b):
    b = mp.mpf(b)
    c = -mp.log(b)
    a = mp.lambertw(c).real / c
    lam = a * mp.log(b)
    # phi(h) = f(a+h) - a = a (e^{-c h} - 1)
    phi = [mp.mpf(0)] + [a * (-c) ** k / mp.factorial(k) for k in range(1, NSER)]
    # Potenzen von phi als Reihen
    def mul(u, v):
        w = [mp.mpf(0)] * NSER
        for i, ui in enumerate(u):
            if ui == 0:
                continue
            for j in range(NSER - i):
                w[i + j] += ui * v[j]
        return w
    pw = [None, phi]
    for k in range(2, NSER):
        pw.append(mul(pw[-1], phi))
    s = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (NSER - 2)
    for n in range(2, NSER):
        acc = mp.mpf(0)
        for k in range(1, n):
            acc += s[k] * pw[k][n]
        s[n] = acc / (lam - lam ** n)
    return dict(b=b, c=c, a=a, lam=lam, s=s, lnb=mp.log(b))


def f(P, z):
    return mp.exp(P['lnb'] * z)


def finv(P, z):
    return mp.log(z) / P['lnb']


def sig_loc(P, h):
    return mp.polyval(P['s'][::-1], h)


def sigma(P, z, tol=mp.mpf('1e-3')):
    n = 0
    while abs(z - P['a']) > tol:
        z = f(P, z)
        n += 1
        if n > 400000:
            raise RuntimeError('kein Einlauf')
    return sig_loc(P, z - P['a']) / P['lam'] ** n


def psi_loc(P, w):
    h = w
    for _ in range(80):
        F = sig_loc(P, h) - w
        dF = mp.polyval([k * P['s'][k] for k in range(len(P['s']) - 1, 0, -1)], h)
        dh = F / dF
        h -= dh
        if abs(dh) < mp.mpf(10) ** (-mp.mp.dps + 3):
            break
    return P['a'] + h


def psi(P, w, tol=mp.mpf('1e-3')):
    n = 0
    while abs(w * P['lam'] ** n) > tol:
        n += 1
    z = psi_loc(P, w * P['lam'] ** n)
    for _ in range(n):
        z = finv(P, z)
    return z


def T_p(P, x, p):
    return psi(P, abs(P['lam']) ** x * p(x))


def T_cos(P, x):
    s1 = sigma(P, mp.mpf(1))
    return T_p(P, x, lambda t: s1 * mp.cospi(t))


def T_reg(P, z, sign=+1):
    """komplexe regulaere Tetration, log lambda = ln|lambda| + sign*i*pi"""
    s1 = sigma(P, mp.mpf(1))
    ell = mp.log(abs(P['lam'])) + sign * 1j * mp.pi
    w = s1 * mp.exp(z * ell)
    # Psi komplex: gleiche Rueckwaertsiteration mit Hauptzweig-Logarithmus
    n = 0
    while abs(w * P['lam'] ** n) > mp.mpf('1e-3'):
        n += 1
    z0 = psi_loc(P, w * P['lam'] ** n)
    for _ in range(n):
        z0 = mp.log(z0) / P['lnb']
    return z0


if __name__ == '__main__':
    mp.nprint = None
    for b in ['0.1', '0.5', '0.9']:
        P = setup(b)
        s1 = sigma(P, mp.mpf(1))
        print(f'b={b}: a={mp.nstr(P["a"], 25)}  lambda={mp.nstr(P["lam"], 25)}  sigma(1)={mp.nstr(s1, 25)}')
        # Pruefung der Funktionalgleichung
        for x in ['0.25', '0.5', '0.7', '1.3', '-0.5', '-1.5']:
            x = mp.mpf(x)
            t = T_cos(P, x)
            t1 = T_cos(P, x + 1)
            print(f'   T_cos({mp.nstr(x,3)}) = {mp.nstr(t, 30)}   |T(x+1)-b^T(x)| = {mp.nstr(abs(t1 - f(P, t)), 3)}')
        print('   T_cos(1/2) - a =', mp.nstr(T_cos(P, mp.mpf(0.5)) - P['a'], 3))
        tr = T_reg(P, mp.mpf(0.5))
        print('   T_reg(1/2) =', mp.nstr(tr, 20), '  Re T_reg(1/2) - a =', mp.nstr(tr.real - P['a'], 10))
