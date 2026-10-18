"""Periodenmittel auf der reellen Achse: Multiplikativitaet und Kennzeichen (Etappe-3-Formeln) fuer den Zoo.

Periode x in [-1, 0] (S(-1) = v, S(0) = f(v)); Substitution x = -1 + sin^2(pi s/2) glaettet die
Wurzel-Singularitaeten an beiden Enden.  Gegenprobe: T(x) = S(x + eps sin 2 pi x).
"""
import numpy as np, json, os, sys, math
from zoo import *
from s1_loesen import ABB

def S_reell(s, x):
    v = s.Sweit(np.asarray(x, complex))
    return v.real + 0j

def Sp_reell(s, x, r=0.02, K=24):
    w = np.exp(2j * np.pi * np.arange(K) / K)
    z = np.asarray(x, float)[:, None] + r * w[None, :]
    v = s.Sweit(z.ravel()).reshape(z.shape)
    return ((v * np.conj(w)[None, :]).mean(1) / r).real

def theta_daten(a, s, x, eps=0.0):
    phi = x + eps * np.sin(2 * np.pi * x); dphi = 1 + 2 * np.pi * eps * np.cos(2 * np.pi * x)
    S = S_reell(s, phi); Sp = Sp_reell(s, phi) * dphi
    ls, n = a.logsigma(S)
    im = np.unwrap(ls.imag)
    h = (ls.real + 1j * im) / a.loglam
    th = h - x
    w = Sp * a.hp(S)          # 1 + theta'
    return th, w, S

res = {}
for nm in (sys.argv[1:] or ['exp', 'q1', 'cosh', 'q4']):
    a = ABB[nm][0]()
    out = {}
    for e in 'AB':
        s = Loesung.laden(a, f'{nm}_{e}')
        for Nq in [200, 300]:
            u, wq = np.polynomial.legendre.leggauss(Nq)
            sq = (u + 1) / 2; wq = wq / 2
            x = -1 + np.sin(np.pi * sq / 2) ** 2; dx = (np.pi / 2) * np.sin(np.pi * sq)
            for eps in [0.0, 1e-4]:
                th, w, S = theta_daten(a, s, x, eps)
                M = lambda G: np.sum(wq * dx * G)
                m1 = M(th)
                lw = np.log(w); lw = lw.real + 1j * np.unwrap(lw.imag)
                lw = lw - 2j * np.pi * np.round(lw[len(lw) // 2].imag / (2 * np.pi))
                d = dict(M_theta=m1,
                         V2=M(th ** 2) - m1 ** 2,
                         V3=M((th - m1) ** 3),
                         expmom=M(np.exp(1j * th)) / np.exp(1j * m1) - 1,
                         sigma_prod=M(np.exp(a.loglam * (th - m1))) * M(np.exp(-a.loglam * (th - m1))) - 1,
                         inv_w=M(1 / w) - 1,
                         log_w=M(lw),
                         int_lnSp=M(np.log(Sp_reell(s, x + eps * np.sin(2 * np.pi * x)) * (1 + 2 * np.pi * eps * np.cos(2 * np.pi * x)))))
                for m in [1, 2, 3]:
                    d[f'neg{m}'] = M(th * np.exp(2j * np.pi * m * x))
                out[(e, Nq, eps)] = d
    print(f'== {nm}')
    keys = list(out[('A', 300, 0.0)].keys())
    for k in keys:
        v0 = out[('A', 300, 0.0)][k]; v1 = out[('B', 200, 0.0)][k]; ve = out[('A', 300, 1e-4)][k]
        print(f'  {k:11s} Kneser: {abs(v0) if k not in ("M_theta","int_lnSp") else v0:.3e}   Gegenprobe A300/B200: {abs(v0-v1):.1e}   eps=1e-4: {abs(ve) if k not in ("M_theta","int_lnSp") else abs(ve-v0):.2e}')
    res[nm] = {k: [complex(v).real, complex(v).imag] for k, v in out[('A', 300, 0.0)].items()}
    res[nm + '_eps'] = {k: [complex(v).real, complex(v).imag] for k, v in out[('A', 300, 1e-4)].items()}
json.dump(res, open(os.path.join(DATEN, 'mittel.json'), 'w'), indent=1)
