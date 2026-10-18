"""Die Boettcher-Superfunktion von f = x^d + c (zweite natuerliche reelle Loesung) gegen Kneser.

phi(f(w)) = phi(w)^d, phi(w) ~ w;  S_B(z) = phi^{-1}(exp(G(c) d^z)), S_B(0) = c.
log phi(w) = log w + sum_n d^{-(n+1)} log(1 + c / f^n(w)^d)
Kneser (Loeser): S_K(-1) = c.  Vergleich: S_B(x) = S_K(x - 1 + p(x)), p 1-periodisch.
Defekt: theta_B(z) = h(S_B(z)) - z hat negative Moden; N = sum |hat(-m)|^2, E = 2 sum m |hat(-m)|^2.
"""
import numpy as np, json, os, sys, math
from zoo import *
from s1_loesen import ABB

def logphi(w, c, d, nmax=60):
    w = np.asarray(w, complex)
    lp = np.log(w); dlp = 1 / w
    fn = w.copy(); dfn = np.ones_like(w)
    for n in range(nmax):
        q = c / fn ** d
        lp = lp + np.log1p(q) / d ** (n + 1)
        dq = -d * q / fn * dfn
        dlp = dlp + dq / (1 + q) / d ** (n + 1)
        if np.all(abs(q) < 1e-18): break
        dfn = d * fn ** (d - 1) * dfn; fn = fn ** d + c
        if np.any(abs(fn) > 1e150): break
    return lp, dlp

def green(c, d):
    lp, _ = logphi(np.array([c + 0j]), c, d)
    return lp[0].real

def SB(z, w0, c, d, Gc):
    z = np.asarray(z, complex); w = np.asarray(w0, complex).copy()
    ziel = Gc * d ** z
    for _ in range(60):
        lp, dlp = logphi(w, c, d)
        dw = (lp - ziel) / dlp; w = w - dw
        if np.max(abs(dw)) < 1e-16 * np.max(abs(w)): break
    return w

def AK(s, w, z0):
    """Kneser-Abelfunktion: Loesung von S_K(z) = w nahe z0."""
    z = np.asarray(z0, complex).copy()
    for _ in range(40):
        v = s.Sweit(z); dv = (s.Sweit(z + 1e-5) - s.Sweit(z - 1e-5)) / 2e-5
        dz = (v - w) / dv; z = z - dz
        if np.max(abs(dz)) < 1e-15: break
    return z

if __name__ == '__main__':
    erg = {}
    for nm, d in [('q1', 2), ('q4', 4)]:
        a = ABB[nm][0](); c = 1.0
        Gc = green(c, d)
        r = dict(G_c=Gc)
        for e in 'AB':
            s = Loesung.laden(a, f'{nm}_{e}')
            M = 512
            x = np.arange(M) / M
            wK = s.Sweit(x - 1 + 0j)
            wB = SB(x, wK, c, d, Gc)
            assert abs(wB[0] - c) < 1e-12, wB[0]
            zK = AK(s, wB.real + 0j, x - 1 + 0j)
            p = (zK.real + 1 - x)
            ph = np.fft.fft(p) / M
            r[e] = dict(maxp=float(abs(p).max()), p_hat=[[v.real, v.imag] for v in ph[:8]])
            # Defekt
            neg = {}
            for Y in [0.03, 0.06]:
                Mz = 4096; xx = np.arange(Mz) / Mz; z = xx + 1j * Y
                w0 = s.Sweit(z - 1 + 0j)   # Startwert: Kneser
                wb = w0.copy()
                # Fortsetzung in Y in kleinen Schritten
                for yy in np.linspace(0, Y, 7)[1:]:
                    wb = SB(xx + 1j * yy, wb, c, d, Gc)
                assert np.all(wb.imag > -1e-14), wb.imag.min()
                ls, _ = a.logsigma(wb)
                im = np.unwrap(np.concatenate([ls.imag, ls.imag[:1] ]))
                th = (ls.real + 1j * im[:-1]) / a.loglam - z
                sp = abs((ls.real[0] + 1j * im[-1]) / a.loglam - (z[0] + 1) - th[0])  # Periodizitaet
                cf = np.fft.fft(th) / Mz
                m = np.arange(1, 40)
                hneg = cf[-m] * np.exp(-2 * np.pi * m * Y)
                neg[Y] = hneg
                r[e][f'per_{Y}'] = float(sp)
            h1, h2 = neg[0.03], neg[0.06]
            r[e]['neg'] = [[v.real, v.imag] for v in h1[:12]]
            r[e]['neg_Ydiff'] = float(np.max(abs(h1[:12] - h2[:12])))
            r[e]['N'] = float(np.sum(abs(h1) ** 2)); r[e]['E'] = float(2 * np.sum(np.arange(1, 40) * abs(h1) ** 2))
            # Vergleich: Mobius-/lineare Naeherung  hat(-m) ~ conj(p_hat(m))*(...) nicht noetig
        print(nm, ' G(c) =', Gc)
        for e in 'AB':
            print(f'  {e}: max|p| = {r[e]["maxp"]:.6e}  p_hat(1) = {complex(*r[e]["p_hat"][1]):.6e}  N = {r[e]["N"]:.6e}  E = {r[e]["E"]:.6e}  '
                  f'|neg(Y=.03)-neg(Y=.06)| = {r[e]["neg_Ydiff"]:.1e}  Per.-Fehler {r[e]["per_0.03"]:.1e}')
        print('   hat theta_B(-m), m=1..6:', ['%.3e' % abs(complex(*v)) for v in r['A']['neg'][:6]])
        erg[nm] = r
    json.dump(erg, open(os.path.join(DATEN, 'boettcher.json'), 'w'), indent=1)
