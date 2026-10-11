"""Strang P — Satz P-V (zertifizierte Kneser-Werte), gemeinsame Bausteine (acb/arb, python-flint).

Bezeichnungen: L = e^L (oberer Fixpunkt), log L = L.  P = Poincaré-Funktion, P(Lw) = e^{P(w)}, P(0) = L, P'(0) = 1.
Obere exakte Lösung Ũ(z) = P(L^{z+θ̃(z)}) = P(exp(L·(z + θ̃(z)))), θ̃(z) = Σ_{k≥0} c_k e^{2πikz}.
Untere Karte T̃(z) = Σ t_n z^n (Taylor von sexp um 0).
Alles hier ist für die NUMERIK (Mittelpunkte); die rigorosen Schranken stehen in eigenen Skripten.
"""
from flint import acb, arb, ctx

_P = {}


def setup(prec, n0=None):
    """Koeffizienten von P bis n0 (Rekursion p_n(L^n − L) = L·[wⁿ]exp(P − L))."""
    ctx.prec = prec
    if n0 is None:
        n0 = int(prec * 0.45) + 40
    L = -acb(-1).lambertw(-1)
    p = [L, acb(1)] + [acb(0)] * (n0 - 1)
    # exp-Reihe inkrementell: E = exp(A), A = P − L; e_m aus n e_n = Σ k a_k e_{n−k}
    a = [acb(0), acb(1)] + [acb(0)] * (n0 - 1)
    e = [acb(1)] + [acb(0)] * n0
    e[1] = acb(1)
    Lpow = L
    for n in range(2, n0 + 1):
        Lpow = Lpow * L
        # e_n mit a_n = 0 vorläufig
        s = acb(0)
        for k in range(1, n):
            s += k * a[k] * e[n - k]
        q = s / n                       # [w^n] E bei a_n = 0
        pn = L * q / (Lpow - L)
        pn = acb(pn.real.mid(), pn.imag.mid())
        p[n] = pn
        a[n] = pn
        e[n] = q + pn                   # a_n trägt mit n·a_n·e_0/n = a_n bei
    _P.update(L=L, p=p, n0=n0, r0=arb("0.12"), prec=prec)
    return L


def L():
    return _P["L"]


def P_dP(w):
    """(P(w), P'(w)) numerisch: w → w/L^k bis |w| ≤ r0, Taylor, dann k-mal exp mit Kettenregel."""
    Lc = _P["L"]
    p = _P["p"]
    k = 0
    while abs(w).mid() > _P["r0"].mid():
        w = w / Lc
        k += 1
    # Horner für P und P'
    s = acb(0)
    ds = acb(0)
    for n in range(len(p) - 1, -1, -1):
        ds = ds * w + s
        s = s * w + p[n]
    # ds ist jetzt P'(w) (Horner-Ableitung)
    for _ in range(k):
        e = s.exp()
        ds = ds * e / Lc               # d/dw [P(Lw)] = L P'(Lw) = e^{P(w)} P'(w)
        s = e
    return s, ds


def Psi(y, w0=None, schritte=60):
    """Ψ = P⁻¹ (Schröder-Funktion) numerisch: logs bis nahe L, dann Newton, Ψ(y) = L^m Ψ(log^m y)."""
    Lc = _P["L"]
    m = 0
    while abs(y - Lc).mid() > arb("0.05").mid():
        y = y.log()
        m += 1
        if m > 400:
            raise RuntimeError("Psi: keine Konvergenz der log-Iteration")
    w = y - Lc if w0 is None else w0
    for _ in range(schritte):
        f, df = P_dP(w)
        dw = (f - y) / df
        w = w - dw
        if abs(dw).mid() < arb(2) ** (-ctx.prec + 20):
            break
    return w, m                          # Ψ(y_orig) = L^m · w


def theta_eval(c, z, mit_abl=False):
    """θ̃(z) = Σ c_k q^k, q = e^{2πiz}; optional θ̃'(z)."""
    q = (2 * acb.pi() * acb(0, 1) * z).exp()
    s, ds = acb(0), acb(0)
    for k in range(len(c) - 1, -1, -1):
        s = s * q + c[k]
        if mit_abl:
            ds = ds * q + k * c[k]
    if not mit_abl:
        return s
    return s, ds * 2 * acb.pi() * acb(0, 1)


def U_dU(c, z):
    """Ũ(z), Ũ'(z) mit Ũ = P(exp(L(z + θ̃(z))))."""
    Lc = _P["L"]
    th, dth = theta_eval(c, z, mit_abl=True)
    w = (Lc * (z + th)).exp()
    f, df = P_dP(w)
    return f, df * w * Lc * (1 + dth)


def T_dT(t, z):
    s, ds = acb(0), acb(0)
    for n in range(len(t) - 1, -1, -1):
        ds = ds * z + s
        s = s * z + t[n]
    return s, ds


def mid(x):
    return acb(x.real.mid(), x.imag.mid())
