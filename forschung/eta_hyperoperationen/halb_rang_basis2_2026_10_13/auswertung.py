"""Auswertung: Halbrang-Jet u_j = c_{1/2}^{(j)}, Rang 3/2, Rangrekursion, Summen bei eps=1.

Rangfortsetzung: Newtonreihe im Rang (Auswahl polynomialen Wachstums), beschleunigt mit
Levin-u.  Fehlerschaetzung: Abstand der Levin-Werte mit N und N-8 Rangen.

Aufruf:  python auswertung.py N40J8R7
"""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from halbschritt import lade, newton_partial

mp.mp.dps = 120


def fortsetzen(seq, s):
    """Wert der Rangfortsetzung bei s (Levin-u auf Newton-Partialsummen)."""
    ps = newton_partial(seq, s)
    terms = [ps[0]] + [b - a for a, b in zip(ps, ps[1:])]
    skala = max(abs(t) for t in terms) or 1
    # endliche Newtonreihe (Polynom im Rang): Abbruch
    if all(abs(t) < skala * mp.mpf(10) ** -80 for t in terms[len(terms) // 2:]):
        return ps[-1]
    while len(ps) > 1 and abs(ps[0]) < skala * mp.mpf(10) ** -90:
        ps = ps[1:]
    L = mp.levin(method='levin', variant='u')
    v, e = L.update_psum(ps)
    return v


def wert(d, p, j, r, s, N=None):
    seq = [mp.mpc(mp.mpf(row[r][j][0]), mp.mpf(row[r][j][1])) for row in d['werte'][p]]
    if N is not None:
        seq = seq[:N + 1]
    return fortsetzen(seq, s)


def wert_ganz(d, k, j, r, s, N=None):
    seq = [mp.mpf(row[r][j]) for row in d['ganz'][str(k)]]
    if N is not None:
        seq = seq[:N + 1]
    return fortsetzen(seq, s)


def reihe_mul(a, b, J):
    c = [mp.mpc(0)] * (J + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j in range(J + 1 - i):
            c[i + j] += x * b[j]
    return c


def rekursions_residuum(Ux1, Wx, Wx1, J):
    """Residuum W(x+1) - U(W(x)) als eps-Reihe.  Ux1[r][j]: r-te Ableitung von u_j bei x+1;
    Wx[j], Wx1[j]: Jet von w=W-(x+1) bei x bzw. x+1 (Index j=1..J, Index 0 = 0)."""
    w = [mp.mpc(0)] + list(Wx)
    pred = list(w)  # U(x+1+w) = x+2+w + sum_r u^{(r)}(x+1) w^r/r!
    pw = [mp.mpc(1)] + [mp.mpc(0)] * J
    for r in range(J):
        for j in range(1, J + 1):
            c = Ux1[r][j - 1] / mp.factorial(r)
            for i in range(J + 1 - j):
                pred[i + j] += c * pw[i]
        pw = reihe_mul(pw, w, J)
    ist = [mp.mpc(0)] + list(Wx1)
    return [ist[j] - pred[j] for j in range(1, J + 1)]


def main(name):
    d = lade(name)
    J = d['J']; R = d['R']; N = d['N']
    out = {'quelle': name, 'J': J, 'N': N}
    # 1. Halbrang- und 3/2-Jets an allen Punkten (Ordnung j, Ableitung r)
    jets = {}
    for p in d['punkte']:
        jets[p] = {}
        for s in ['0.5', '1.5']:
            tab = []; err = []
            for r in range(R + 1):
                row = []; erow = []
                for j in range(J):
                    v = wert(d, p, j, r, mp.mpf(s))
                    v2 = wert(d, p, j, r, mp.mpf(s), N=N - 8)
                    row.append(v); erow.append(abs(v - v2))
                tab.append(row); err.append(erow)
            jets[p][s] = (tab, err)
        tab, err = jets[p]['0.5']
        print('x =', p, ' u_j (r=0):', ' '.join('%s(+-%s)' % (mp.nstr(mp.re(tab[0][j]), 6), mp.nstr(err[0][j], 1)) for j in range(J)))
    out['halbrang_jet'] = {p: {'u': [[str(mp.re(v)), str(mp.im(v))] for v in jets[p]['0.5'][0][0]],
                               'fehler': [float(e) for e in jets[p]['0.5'][1][0]],
                               'u32': [[str(mp.re(v)), str(mp.im(v))] for v in jets[p]['1.5'][0][0]],
                               'fehler32': [float(e) for e in jets[p]['1.5'][1][0]]} for p in d['punkte']}
    # 2. Rangrekursion W(x+1)=U(W(x)) mit U=T^{1/2}f, W=T^{3/2}f, Kontrolle mit ganzzahligen Raengen
    paare = [('0.25', '1.25'), ('0.5', '1.5'), ('0.75', '1.75'), ('0.5+0.5j', '1.5+0.5j'), ('0.25+0.25j', '1.25+0.25j'),
             ('1.5', '2.5'), ('1.25', '2.25'), ('-0.5', '0.5')]
    rek = {}
    for x, x1 in paare:
        Ux1 = jets[x1]['0.5'][0]
        Wx = jets[x]['1.5'][0][0]; Wx1 = jets[x1]['1.5'][0][0]
        res = rekursions_residuum(Ux1, Wx, Wx1, J)
        # Kontrolle: ganzzahlige Raenge n=1 (U=T^0 f=f, W=T^1 f) direkt aus den Daten
        Uc = [[mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in d['werte'][x1][0][r]] for r in range(R + 1)]
        Wc = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in d['werte'][x][1][0]]
        Wc1 = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in d['werte'][x1][1][0]]
        res_int = rekursions_residuum(Uc, Wc, Wc1, J)
        skala = [abs(v) for v in Wx1]
        rek[x] = {'residuum': [float(abs(v)) for v in res], 'betrag_W': [float(v) for v in skala],
                  'kontrolle_ganzzahlig': [float(abs(v)) for v in res_int]}
        print('Rangrekursion x=%s:' % x, ' '.join(mp.nstr(abs(v), 2) for v in res), ' | ganzzahlig:', mp.nstr(max(abs(v) for v in res_int), 2))
    out['rangrekursion'] = rek
    # 3. Ganzzahlige Hoehen: T(U f)(k) gegen T^{3/2} f(k)
    K = d['K']
    U = {k: [[wert_ganz(d, k, j, r, mp.mpf('0.5')) for j in range(J)] for r in range(J)] for k in range(K + 1)}
    W32 = {k: [wert_ganz(d, k, j, 0, mp.mpf('1.5')) for j in range(J)] for k in range(K + 1)}
    ganz = {}
    e = [mp.mpc(0)] * J  # e(0)=0 fuer S=T(Uf)
    for k in range(K):
        # S(k+1) = U(S(k)) = U(k+1+e(k))
        res_e = rekursions_residuum([[mp.mpc(v) for v in row] for row in U[k + 1]], e, [mp.mpc(0)] * J, J)
        e = [-v for v in res_e]  # ist=0 -> res = -pred, also pred = -res
        diff = [abs(e[j] - W32[k + 1][j]) for j in range(J)]
        ganz[k + 1] = {'T_von_U': [str(mp.re(v)) for v in e], 'T32': [str(v) for v in W32[k + 1]], 'abstand': [float(v) for v in diff]}
        print('ganze Hoehe %d: |T(Uf)-T^{3/2}f| je Ordnung:' % (k + 1), ' '.join(mp.nstr(v, 2) for v in diff))
    out['ganze_hoehen'] = ganz
    print('Halbrang an ganzen Hoehen (Werte u_j(k)):')
    for k in range(K + 1):
        print('  k=%d' % k, ' '.join(mp.nstr(U[k][0][j], 8) for j in range(J)))
    (HIER / 'data' / f'auswertung_{name}.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'N40J8R7')
