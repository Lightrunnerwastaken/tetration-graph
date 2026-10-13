"""Exponentialpolynome  sum_m 2^(m x) P_m(x)  mit Ballarithmetik (python-flint arb_poly).

Grundbausteine fuer den Rangschritt T(f)=S, S(0)=1, S(x+1)=f(S(x)) als Jet in eps um den
Successor x+1.  Auswahl: die Summation (Antidifferenz) wird in der Klasse der Exponential-
polynome geloest; das ist die eindeutige Loesung mit polynomialem vertikalem Wachstum
(siehe MATHEMATIK.md, Paragraph 2).
"""
from flint import arb_poly, arb, acb, ctx, fmpq_poly, fmpq

LN2 = None
SHIFT1 = None
_bern = {}


def setprec(bits):
    global LN2, SHIFT1
    ctx.prec = bits
    LN2 = arb(2).log()
    SHIFT1 = arb_poly([1, 1])
    _bern.clear()


setprec(400)


def _sumpoly(k):
    """Polynom q_k mit q_k(x+1)-q_k(x)=x^k, q_k(0)=0:  (B_{k+1}(x)-B_{k+1})/(k+1)."""
    if k not in _bern:
        b = fmpq_poly.bernoulli_poly(k + 1)
        c = [fmpq(q) for q in b.coeffs()]
        c[0] = fmpq(0)
        _bern[k] = arb_poly([arb(int(q.p)) / arb(int(q.q)) / (k + 1) for q in c])
    return _bern[k]


def clean(e):
    return {m: p for m, p in e.items() if p.length() > 0}


def add(a, b, s=1):
    out = dict(a)
    for m, p in b.items():
        out[m] = out[m] + p * s if m in out else p * s
    return clean(out)


def scale(a, c):
    return clean({m: p * c for m, p in a.items()})


def mul(a, b):
    out = {}
    for m1, p1 in a.items():
        for m2, p2 in b.items():
            m = m1 + m2
            out[m] = out[m] + p1 * p2 if m in out else p1 * p2
    return clean(out)


def shift(a, h=1):
    """f(x+h) fuer ganzzahliges h>=0 (hier nur h=1 benutzt)."""
    assert h == 1
    return clean({m: (p(SHIFT1)) * (arb(2) ** m) for m, p in a.items()})


def deriv(a):
    out = {}
    for m, p in a.items():
        q = p.derivative()
        if m:
            q = q + p * (m * LN2)
        out[m] = q
    return clean(out)


def antidiff(a):
    """Loesung e von e(x+1)-e(x)=a(x) mit e(0)=0, im Exponentialpolynom-Raum."""
    out = {}
    for m, p in a.items():
        if m == 0:
            q = arb_poly([0])
            for k, c in enumerate(p.coeffs()):
                q = q + _sumpoly(k) * c
            out[0] = q
        else:
            # lam*q(x+1)-q(x)=p(x): Rueckwaertssubstitution in der Monombasis
            lam = arb(2) ** m
            pc = p.coeffs()
            D = len(pc) - 1
            c = [arb(0)] * (D + 1)
            inv = 1 / (lam - 1)
            for k in range(D, -1, -1):
                acc = pc[k]
                b = 1
                for i in range(k + 1, D + 1):
                    b = b * i // (i - k)  # binom(i,k)
                    acc -= lam * c[i] * b
                c[k] = acc * inv
            q = arb_poly(c)
            out[m] = q
    out = clean(out)
    c0 = evaluate(out, 0)
    return add(out, {0: arb_poly([-c0.real])})


def evaluate(a, x):
    """Wert bei reellem oder komplexem x (acb)."""
    x = acb(x)
    tot = acb(0)
    for m, p in a.items():
        tot += p(x) * (acb(2) ** (m * x) if m else 1)
    return tot


def degree(a):
    return max([p.degree() for p in a.values()] + [-1])


def d_base():
    """d(x)=2^x-x-1."""
    return {0: arb_poly([-1, -1]), 1: arb_poly([1])}


def rank_step(g, J):
    """Ein voller Rangschritt T auf dem Jet g=[g_1..g_J] (Ordnung eps^j), gibt e=[e_1..e_J].

    e(x+1)-e(x) = delta(x+1+e(x)),  delta=sum eps^i g_i,  e(0)=0.
    """
    # Ableitungen der verschobenen Eingabe: gd[i][r] = g_i^{(r)}(x+1)/r!
    gd = []
    for i in range(J):
        row = []
        h = shift(g[i]) if g[i] else {}
        fact = 1
        for r in range(J - i):
            row.append(scale(h, arb(1) / fact) if h else {})
            h = deriv(h) if h else {}
            fact *= (r + 1)
        gd.append(row)
    e = [dict() for _ in range(J)]
    # pw[r][k] = [eps^k] e^r, inkrementell: Eintrag k wird bei Schritt j=k+1 gebildet
    pw = [[{0: arb_poly([1])}] + [{} for _ in range(J)]] + [[{} for _ in range(J + 1)] for _ in range(J)]
    for j in range(1, J + 1):
        k = j - 1
        for r in range(1, k + 1):
            acc = {}
            for l in range(1, k - r + 2):
                if e[l - 1] and pw[r - 1][k - l]:
                    acc = add(acc, mul(e[l - 1], pw[r - 1][k - l]))
            pw[r][k] = acc
        G = {}
        for i in range(1, j + 1):
            kk = j - i
            for r in range(0, kk + 1):
                term = gd[i - 1][r] if r < len(gd[i - 1]) else {}
                if term and pw[r][kk]:
                    G = add(G, mul(term, pw[r][kk]))
        e[j - 1] = antidiff(G) if G else {}
    return e
