"""Kanonische Rangleiter fuer kleine Basen 1 < b < e^(1/e).

Fuer b < e^(1/e) hat jede Stufe H_r(b,.) (r >= 3) einen anziehenden Fixpunkt p_r >= 1
(Satz M1/M2 im Graphen), und die naechste Stufe ist die Koenigs-Superfunktion
    H_{r+1}(x) = sigma_r^{-1}( lambda_r^x * sigma_r(1) ),   lambda_r = H_r'(p_r),
mit der Schroeder-Funktion sigma_r von H_r am Fixpunkt p_r. Das ist die einzige
reell-analytische, monotone, beschraenkte Loesung (Koenigs-Eindeutigkeit fuer
x -> +inf gegen p_r mit Rate lambda^x ohne periodischen Anteil).

Jede Stufe wird auf [0, X] durch ein Tschebyschow-Polynom (mpmath, dps einstellbar)
dargestellt. Aufruf:
    python koenigs_leiter.py b R dps grad   -> json in daten/
"""
import sys, json, os
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


class Cheb:
    def __init__(self, f, a, b, n):
        self.a, self.b, self.n = mp.mpf(a), mp.mpf(b), n
        xs = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / n) for k in range(n)]
        fx = [f(self.a + (self.b - self.a) * (t + 1) / 2) for t in xs]
        self.c = []
        for j in range(n):
            s = mp.fsum(fx[k] * mp.cos(mp.pi * j * (k + mp.mpf(1) / 2) / n) for k in range(n))
            self.c.append(2 * s / n)
        self.c[0] /= 2

    def __call__(self, x):
        t = (2 * x - self.a - self.b) / (self.b - self.a)
        b1 = b2 = mp.mpf(0)
        for cj in reversed(self.c[1:]):
            b1, b2 = 2 * t * b1 - b2 + cj, b1
        return t * b1 - b2 + self.c[0]

    def deriv(self, x, k=1):
        return mp.diff(self, x, k)


def leiter(b, R, X=2, grad=40, tol=None, eps_theta=0, seed4=None, X0=0):
    b = mp.mpf(b)
    tol = tol or mp.mpf(10) ** (-(mp.mp.dps - 5))
    H = {3: (lambda x: b ** x)}
    start = 3
    if seed4 is not None:
        H[4] = seed4
        start = 4
    out = []
    for r in range(start, R):
        f = H[r]
        # Fixpunkt p_r: Iteration von 1 aus, dann Newton
        y = mp.mpf(1)
        for _ in range(4000):
            y2 = f(y)
            if abs(y2 - y) < mp.mpf(10) ** (-12):
                break
            y = y2
        for _ in range(60):
            g = f(y) - y
            d = mp.diff(f, y) - 1
            y2 = y - g / d
            if abs(y2 - y) < tol:
                y = y2
                break
            y = y2
        p = y
        lam = mp.diff(f, p)
        f2 = mp.diff(f, p, 2)
        f3 = mp.diff(f, p, 3)
        # sigma(p+u) = u + c2 u^2 + c3 u^3 + ...
        c2 = f2 / (2 * lam * (1 - lam))
        # c3 aus sigma(f(p+u)) = lam sigma(p+u), Ordnung u^3
        c3 = (f3 / 6 + 2 * c2 * lam * f2 / 2) / (lam - lam ** 3)
        eps = mp.mpf(10) ** (-(mp.mp.dps / 4.0))

        def sigma(yv):
            n = 0
            u = yv - p
            z = yv
            while abs(u) > eps:
                z = f(z)
                u = z - p
                n += 1
                if n > 100000:
                    raise RuntimeError('keine Konvergenz')
            return (u + c2 * u ** 2 + c3 * u ** 3) / lam ** n

        s1 = sigma(mp.mpf(1))

        def Hn(x):
            w = lam ** x * s1
            # sigma^{-1}(w): Start lineare Naeherung ueber Rueckwaertsiteration
            # loese sigma(y) = w mit Sekante/Bisektion auf [1, p)
            if abs(w) < eps:
                u = w
                for _ in range(8):
                    u = w - c2 * u ** 2 - c3 * u ** 3
                return p + u
            lo, hi = mp.mpf(1) - mp.mpf(1) * (x < 0), p
            # Bisektion grob, dann Sekante
            return mp.findroot(lambda yy: sigma(yy) - w, (lo + (hi - lo) * 0.3, lo + (hi - lo) * 0.7), solver='anderson')

        ch = Cheb(Hn, X0, X, grad)
        if r + 1 == 4 and eps_theta:
            ch0 = ch
            if str(eps_theta).startswith('c'):   # Kosinus-Richtung psi = x + e(1-cos 2 pi x)
                et = mp.mpf(str(eps_theta)[1:])
                ch = Cheb(lambda x: ch0(x + et * (1 - mp.cos(2 * mp.pi * x))), X0, X, grad)
            else:
                et = mp.mpf(eps_theta)
                ch = Cheb(lambda x: ch0(x + et * mp.sin(2 * mp.pi * x)), X0, X, grad)
        H[r + 1] = ch
        # Kontrollen
        res = max(abs(ch(mp.mpf(k) / 8 + 1) - f(ch(mp.mpf(k) / 8))) for k in range(0, 9))
        row = dict(r=r, p=str(p), lam=str(lam), f2=str(f2),
                   next_half=str(ch(mp.mpf('0.5'))), next_1p5=str(ch(mp.mpf('1.5'))),
                   next_2=str(ch(mp.mpf(2))), next_quarter=str(ch(mp.mpf('0.25'))),
                   next_3q=str(ch(mp.mpf('0.75'))), resid=mp.nstr(res, 3))
        if X0 < 0:
            row.update({f'neg_{t}': str(ch(mp.mpf(t))) for t in ('-0.75', '-0.5', '-0.25')})
            row['resid_neg'] = mp.nstr(max(abs(ch(mp.mpf(kk) / 8) - f(ch(mp.mpf(kk) / 8 - 1))) for kk in range(1, 8)), 3)
        out.append(row)
        print(r, mp.nstr(p, 15), mp.nstr(lam, 12), mp.nstr(ch(mp.mpf('0.5')), 15),
              mp.nstr(ch(mp.mpf(2)), 15), 'res', mp.nstr(res, 3), flush=True)
    return out, H


if __name__ == '__main__':
    b = sys.argv[1]
    R = int(sys.argv[2])
    mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    grad = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    et = sys.argv[5] if len(sys.argv) > 5 else '0'
    out, _ = leiter(mp.mpf(eval(b, {'sqrt': mp.sqrt, 'e': mp.e, 'exp': mp.exp})), R, grad=grad, eps_theta=et)
    os.makedirs(os.path.join(HERE, 'daten'), exist_ok=True)
    fn = os.path.join(HERE, 'daten', f'koenigs_b{b.replace("/", "_").replace("(", "").replace(")", "")}_R{R}_d{mp.mp.dps}_g{grad}' + (f'_th{et}' if et != '0' else '') + '.json')
    json.dump(dict(b=b, dps=mp.mp.dps, grad=grad, rows=out), open(fn, 'w'), indent=1)
    print('->', fn)
