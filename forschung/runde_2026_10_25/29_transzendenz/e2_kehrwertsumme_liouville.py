# e2: C_b = sum_k 1/(b^^k) fuer ganze b >= 2 — expliziter Kettenbruch per Faltungslemma, Liouville-Eigenschaft.
# Regel (zu pruefen): S_N = sum_{k<=N} b^{-e_k}, e_1 = 1, e_{k+1} = b^{e_k} (also b^^k = b^{e_k}), reduzierter Nenner q_N = b^{e_N}.
# Ist CF(S_N) = [0; a_1..a_n] mit n gerade (sonst a_n -> a_n - 1, 1), x = b^{e_{N+1} - 2 e_N}, dann
#   CF(S_{N+1}) = [0; a_1..a_n, x - 1, 1, a_n - 1, a_{n-1}, ..., a_1]   (Nullen werden verschmolzen).
from fractions import Fraction
import math, sys

def cf(fr):
    out = []
    p, q = fr.numerator, fr.denominator
    while q:
        a = p // q
        out.append(a)
        p, q = q, p - a * q
    return out

def normalize(a):
    # verschmilzt [.., u, 0, v, ..] -> [.., u+v, ..]; entfernt abschliessende Probleme
    a = list(a)
    i = 1
    while i < len(a):
        if a[i] == 0 and 0 < i < len(a) - 1:
            a[i - 1:i + 2] = [a[i - 1] + a[i + 1]]
            i = max(1, i - 1)
        else:
            i += 1
    if len(a) > 1 and a[-1] == 1:      # kanonische Form: letzter Teilnenner > 1
        a[-2:] = [a[-2] + 1]
    return a

def fold(a, x):
    a = list(a)
    if (len(a) - 1) % 2 == 1:          # n ungerade -> gleichwertige Form mit gerader Laenge
        a[-1:] = [a[-1] - 1, 1]
    tail = [x - 1, 1, a[-1] - 1] + a[-2:0:-1]
    return normalize(a + tail)

def evalcf(a):
    v = Fraction(a[-1])
    for t in reversed(a[:-1]):
        v = t + 1 / v
    return v

def run(b, N):
    e = [1]
    while len(e) < N:
        e.append(b ** e[-1])
    S = Fraction(0)
    pred = None
    ok = True
    for k in range(N):
        S += Fraction(1, b ** e[k])
        c = cf(S)
        assert S.denominator == b ** e[k]
        if pred is not None:
            ok &= (pred == c)
        if k + 1 < N:
            pred = fold(c, b ** (e[k + 1] - 2 * e[k]))
    return S, cf(S), ok, e

if __name__ == '__main__':
    for b, N in [(2, 5), (3, 3), (4, 3), (5, 3), (7, 3), (10, 2)]:
        S, c, ok, e = run(b, N)
        short = [t if t < 10**12 else 'b^%d' % round(math.log(t, b)) if abs(math.log(t, b) - round(math.log(t, b))) < 1e-9 else '~10^%d' % int(math.log10(t)) for t in c]
        print('b=%d N=%d  Faltungsregel stimmt: %s  Laenge %d' % (b, N, ok, len(c)))
        print('   CF(S_N) =', short)
    # Liouville-Guete der Naeherung p_N/q_N fuer b=2: -log|C - S_N| / log q_N >= (e_{N+1} - 1)/e_N
    e = [1, 2, 4, 16, 65536]
    for N in range(1, 5):
        print('b=2 N=%d: Exponent >= %.1f' % (N, (e[N] - 1) / e[N - 1]))
