# e1: Kettenbruch-Statistik einer ganzen Population von Kneser-Werten (Runde 25.10.2026, Gebiet 29).
# Daten: cache/te_values_520.json, cache/te2_values_520.json (Kopien aus Tetration/Research/arith-half, dps 520, >= 494 wahre Stellen).
# Strenge: Kettenbruch des Intervalls [x - 10^-D, x + 10^-D] (gemeinsamer Anfang beider Endpunkte), D = 470.
# Vergleich: dieselbe Pipeline auf Zufallszahlen (Monte Carlo) und auf klassischen Konstanten.
import json, math, random, sys
from fractions import Fraction
import mpmath as mp

D = 470
mp.mp.dps = 520

def cf_interval(s):
    """s: Dezimalstring. Liefert Teilnenner, die fuer alle Zahlen in [x-eps, x+eps] gleich sind."""
    x = Fraction(s.strip())
    eps = Fraction(1, 10**D)
    lo, hi = x - eps, x + eps
    out = []
    while True:
        a, b = math.floor(lo), math.floor(hi)
        if a != b:
            break
        out.append(a)
        lo, hi = lo - a, hi - a
        if lo == 0 or hi == 0:
            break
        lo, hi = 1 / hi, 1 / lo
    return out

def stats(a):
    """a: Teilnenner ab a0. Liefert q_n, Levy, Chintschin, eff. Irrationalitaetsexponent."""
    q_prev, q = 1, a[1]
    logq = [0.0, math.log(q)]
    for k in a[2:]:
        q_prev, q = q, k * q + q_prev
        logq.append(math.log(q))
    n = len(a) - 1           # Zahl der Teilnenner a1..an
    pq = a[1:]
    khin = math.exp(sum(math.log(t) for t in pq) / n)
    levy = logq[-1] / n
    # mu_i = 1 + log q_{i+1} / log q_i  (|x - p_i/q_i| ~ 1/(q_i q_{i+1}))
    mus = [1 + logq[i + 1] / logq[i] for i in range(20, len(logq) - 1)]
    return dict(n=n, khin=khin, levy=levy, mu_max=max(mus), a_max=max(pq[1:]) if len(pq) > 1 else pq[0],
                rel_max=max(math.log(pq[i]) / logq[i] for i in range(20, len(logq) - 1)))

GK_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 6), (7, 10), (11, 20), (21, 50), (51, 10**9)]
def gk_p(lo, hi):
    # P(lo <= a <= hi) = log2((1+1/lo)/(1+1/(hi+1))) fuer Gauss-Kuzmin
    return math.log2((1 + 1 / lo) / (1 + 1 / (hi + 1)))

def chi2(counts, N):
    return sum((c - N * gk_p(lo, hi)) ** 2 / (N * gk_p(lo, hi)) for c, (lo, hi) in zip(counts, GK_BINS))

def bin_counts(pq):
    c = [0] * len(GK_BINS)
    for t in pq:
        for j, (lo, hi) in enumerate(GK_BINS):
            if lo <= t <= hi:
                c[j] += 1
                break
    return c

def main():
    d = json.load(open('cache/te_values_520.json'))
    d2 = json.load(open('cache/te2_values_520.json'))
    dup = {'0.1', '0.2', '0.4', '0.6', '0.8', '0.9'}   # = 1/10, 1/5, 2/5, 3/5, 4/5, 9/10
    vals = {}
    for k, v in d['sexp_e'].items():
        if k not in dup:
            vals['T_e(%s)' % k] = v
    for k, v in d['slog_e'].items():
        vals['slog_e(%s)' % k] = v
    for k, v in d2['sexp_2'].items():
        vals['T_2(%s)' % k] = v
    for k, v in d2['slog_2'].items():
        vals['slog_2(%s)' % k] = v
    # Gegenprobe der Stellen: Funktionalgleichung zwischen gespeicherten Werten
    chk = []
    for x, y in [('-1/2', '1/2'), ('-1/4', '3/4'), ('-1/3', '2/3'), ('-3/4', '1/4'), ('-2/3', '1/3')]:
        err = abs(mp.exp(mp.mpf(d['sexp_e'][x])) - mp.mpf(d['sexp_e'][y]))
        chk.append((x, y, float(mp.log10(err)) if err > 0 else -999))
    v2 = d2['sexp_2']
    err = abs(mp.power(2, mp.mpf(v2['-1/2'])) - mp.mpf(v2['1/2']))
    chk.append(('2:-1/2', '2:1/2', float(mp.log10(err)) if err > 0 else -999))
    print('Gegenprobe log10|exp(T(x))-T(x+1)|:', chk)

    rows, pooled = [], []
    for name, s in vals.items():
        a = cf_interval(s)
        st = stats(a)
        st['name'] = name
        rows.append(st)
        pooled += a[1:]
    N = len(pooled)
    cnt = bin_counts(pooled)
    c2 = chi2(cnt, N)
    khin_pool = math.exp(sum(math.log(t) for t in pooled) / N)
    # Kontrollen: Zufallszahlen durch dieselbe Pipeline
    random.seed(20261025)
    ctrl = []
    for _ in range(int(sys.argv[1]) if len(sys.argv) > 1 else 200):
        s = '1.' + ''.join(random.choice('0123456789') for _ in range(515))
        a = cf_interval(s)
        ctrl.append(stats(a))
    # Kontroll-Population gleicher Groesse fuer chi2 / Chintschin
    ctrl_pools = []
    for rep in range(40):
        pool = []
        for _ in range(len(vals)):
            s = '1.' + ''.join(random.choice('0123456789') for _ in range(515))
            pool += cf_interval(s)[1:]
        ctrl_pools.append((chi2(bin_counts(pool), len(pool)), math.exp(sum(math.log(t) for t in pool) / len(pool)), len(pool)))
    # klassische Konstanten
    classic = {'pi': mp.pi, 'e': mp.e, 'sqrt2': mp.sqrt(2), 'ln2': mp.log(2), 'zeta3': mp.zeta(3), 'gamma': mp.euler,
               'e^e': mp.e ** mp.e, 'L.re': mp.re(-mp.lambertw(-1))}
    crow = {}
    for k, v in classic.items():
        a = cf_interval(mp.nstr(v, 515, strip_zeros=False))
        st = stats(a)
        st['chi2'] = chi2(bin_counts(a[1:]), len(a) - 1)
        crow[k] = st
    out = dict(D=D, n_werte=len(vals), N_pooled=N, gk_counts=cnt, gk_expected=[N * gk_p(lo, hi) for lo, hi in GK_BINS],
               chi2=c2, chi2_df=len(GK_BINS) - 1, khinchin_pooled=khin_pool, rows=rows,
               ctrl=dict(mu_max=sorted(r['mu_max'] for r in ctrl), rel_max=sorted(r['rel_max'] for r in ctrl),
                         levy=[r['levy'] for r in ctrl], khin=[r['khin'] for r in ctrl], n=[r['n'] for r in ctrl]),
               ctrl_pools=ctrl_pools, classic=crow, gegenprobe=chk)
    json.dump(out, open('cache/e1_ergebnis.json', 'w'), indent=1)
    # Kurzausgabe
    print('Werte:', len(vals), ' Teilnenner gesamt:', N)
    print('GK-Zaehlung:', cnt)
    print('GK-erwartet:', [round(N * gk_p(lo, hi), 1) for lo, hi in GK_BINS])
    print('chi2 = %.2f (df %d)' % (c2, len(GK_BINS) - 1), ' Kontroll-chi2 (40 Pools):', sorted(round(p[0], 1) for p in ctrl_pools))
    print('Chintschin gepoolt = %.4f (K0 = 2.685452); Kontrollen:' % khin_pool, sorted(round(p[1], 4) for p in ctrl_pools))
    import statistics as S
    print('Levy: Kneser mean %.4f sd %.4f | Kontrolle mean %.4f sd %.4f | Soll 1.18657' % (
        S.mean(r['levy'] for r in rows), S.stdev(r['levy'] for r in rows), S.mean(out['ctrl']['levy']), S.stdev(out['ctrl']['levy'])))
    print('n (Teilnenner je Wert): Kneser mean %.1f | Kontrolle mean %.1f' % (S.mean(r['n'] for r in rows), S.mean(out['ctrl']['n'])))
    mm = [r['mu_max'] for r in rows]
    cm = out['ctrl']['mu_max']
    print('mu_max: Kneser max %.4f median %.4f | Kontrolle median %.4f, 95%%-Quantil %.4f, max %.4f' % (
        max(mm), S.median(mm), S.median(cm), cm[int(0.95 * len(cm))], max(cm)))
    for r in sorted(rows, key=lambda r: -r['mu_max'])[:5]:
        print('  ', r['name'], 'n=%d mu_max=%.4f a_max=%d levy=%.4f khin=%.3f' % (r['n'], r['mu_max'], r['a_max'], r['levy'], r['khin']))
    for k, st in crow.items():
        print('  Kontrolle', k, 'n=%d chi2=%.1f khin=%.3f levy=%.4f mu_max=%.4f a_max=%d' % (st['n'], st['chi2'], st['khin'], st['levy'], st['mu_max'], st['a_max']))

if __name__ == '__main__':
    main()
