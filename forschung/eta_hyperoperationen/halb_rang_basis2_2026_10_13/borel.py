"""Borel-Pade-Summation einer eps-Reihe sum c_j eps^j bei eps=1 (bzw. eps)."""
import mpmath as mp


def borel_pade(c, L, M, eps=1, richtung=0.0, info=False):
    """c[0..L+M]: Koeffizienten. Borel-Transformierte B(t)=sum c_j t^j/j!, Pade [L/M],
    Laplace-Integral int_0^inf e^{-t} B(eps t) dt entlang Strahl arg t = richtung."""
    b = [c[j] / mp.factorial(j) for j in range(L + M + 1)]
    p, q = mp.pade(b, L, M)
    P = lambda t: mp.polyval(p[::-1], t) / mp.polyval(q[::-1], t)
    w = mp.exp(1j * richtung)
    f = lambda t: mp.exp(-t * w) * P(eps * t * w) * w
    val = mp.quad(f, [0, 2, 8, 30, mp.inf])
    if info:
        roots = mp.polyroots(q[::-1], maxsteps=300, extraprec=300)
        return val, sorted(roots, key=abs)
    return val
