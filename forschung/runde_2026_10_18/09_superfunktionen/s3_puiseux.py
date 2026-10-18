"""Puiseux-Vorhersage der theta-Koeffizienten mit exakten Binomialkoeffizienten.

theta(z) = 1 - z + h(g(S(z))) nahe z1 = -1 (S(z1) = v kritischer Wert, Grad d).
X = g(S) - w_c,  S - v = A zeta + B zeta^2 + ...,  zeta = z - z1 = log(q)/(2 pi i).
  d = 2 :  a_n = E1 [q^n] w^{1/2} + E3 [q^n] w^{3/2} + O(n^{-7/2}),  w = 1 - q
  d = 4 (x^4+c): a_n = sum_{k=1..3} h_k (A i/2pi)^{k/4} [q^n] w^{k/4} + O(n^{-9/4})
"""
import numpy as np, json, sys, math, os
from scipy.special import gammaln, gamma
from zoo import *
from s1_loesen import ABB

def binkoeff(alpha, n):
    """[q^n](1-q)^alpha = Gamma(n-alpha)/(Gamma(-alpha) Gamma(n+1))  (alpha nicht ganz)."""
    n = np.asarray(n, float)
    sgn = np.sign(gamma(-alpha))
    return sgn * np.exp(gammaln(n - alpha) - gammaln(-alpha) - gammaln(n + 1))

def taylor_kreis(F, z0, r, K=32, m=6):
    w = np.exp(2j * np.pi * np.arange(K) / K)
    v = F(z0 + r * w)
    c = np.fft.fft(v) / K
    return c[:m] / r ** np.arange(m)

GOBEN = {'q1': lambda w: 1j * np.sqrt(1.0 - w + 0j),           # = sqrt(w-c) von oben, analytisch bei 0
         'cosh': lambda w: 1j * np.arccos(w + 0j),
         'q4': lambda w: np.exp(1j * np.pi / 4) * (1.0 - w + 0j) ** 0.25}
def h_glatt(a, w, nm):
    # h(w) = 1 + h(g(w)); g von oben analytisch ueber die reelle Achse fortgesetzt
    ls, n = a.logsigma(GOBEN[nm](w)); ls = ls + a.loglam
    im = np.unwrap(np.concatenate([ls.imag, ls.imag[:1]]))[:-1]
    return (ls.real + 1j * im) / a.loglam

erg = {}
for nm in ['q1', 'cosh', 'q4']:
    a = ABB[nm][0](); s = Loesung.laden(a, f'{nm}_A')
    kr = a.krit; d = kr['d']; aa = kr['a']
    St = taylor_kreis(s.Sweit, -1.0, 0.05)   # S(z1+zeta) = v + A zeta + B zeta^2
    A, B = St[1].real, St[2].real
    hk = taylor_kreis(lambda w: h_glatt(a, w, nm), kr['wc'], 0.05, K=64, m=6)
    th = json.load(open(os.path.join(DATEN, f'theta_{nm}.json')))
    an = np.array([complex(*v) for v in th['a']])
    n = np.arange(len(an)); I = 1j / (2 * np.pi)
    if d == 2:
        delta = (A / 12.0) if nm == 'cosh' else 0.0
        E1 = hk[1] * np.sqrt(A / aa * I)
        E3 = hk[1] * np.sqrt(A / aa * I) * (0.25 + (B / (2 * A) - delta) * I) + hk[3] * (A / aa * I) ** 1.5
        P1 = E1 * binkoeff(0.5, n)
        P2 = P1 + E3 * binkoeff(1.5, n)
        terme = {'E1': E1, 'E3': E3}
    else:
        P1 = hk[1] * (A * I) ** 0.25 * binkoeff(0.25, n)
        P2 = P1 + hk[2] * (A * I) ** 0.5 * binkoeff(0.5, n) + hk[3] * (A * I) ** 0.75 * binkoeff(0.75, n)
        terme = {'h1': hk[1], 'h2': hk[2], 'h3': hk[3]}
    print(f'{nm}: A = {A:.15f}  B = {B:.15f}  h_k = {hk[1]:.12f}, {hk[2]:.12f}, {hk[3]:.12f}')
    print('   n     |a-P1|/|a|   |a-P2|/|a|   n^p |a-P2|/|a|')
    p = 2 if d == 2 else 1
    rows = []
    for k in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        r1 = abs(an[k] - P1[k]) / abs(an[k]); r2 = abs(an[k] - P2[k]) / abs(an[k])
        print(f'{k:5d}  {r1:.3e}   {r2:.3e}   {k**p * r2:.4f}')
        rows.append([k, r1, r2])
    erg[nm] = dict(A=A, B=B, h=[[x.real, x.imag] for x in hk[1:4]], terme={k: [v.real, v.imag] for k, v in terme.items()}, rows=rows)
json.dump(erg, open(os.path.join(DATEN, 'puiseux.json'), 'w'), indent=1)
