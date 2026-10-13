"""Ganzzahlige Raenge n=0..N des Jets T^n(f_eps), f_eps=x+1+eps(2^x-x-1), Ordnungen 1..J.

Speichert fuer jede Testhoehe x und jede Ableitungsordnung r<=R die Werte
  d^r/dx^r c_n^{(j)}(x)   (n=0..N, j=1..J)
als Dezimalzeichenketten (Mittelpunkt) und den groessten Kugelradius (Rundungskontrolle).
Exakt im Raum der Exponentialpolynome; nur Rundung durch Ballarithmetik.

Beispiel:  python raenge.py --N 60 --J 8 --prec 800 --name N60J8
"""
import sys, json, time, argparse
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import exppoly as E
from flint import acb, arb

PUNKTE = ['0.25', '0.5', '0.75', '1.25', '1.5', '1.75', '-0.5', '0.5+0.5j', '1.5+0.5j', '0.25+0.25j', '1.25+0.25j', '0.5+1j', '2.5', '2.25']


def lograd(v):
    r = v.rad()
    if r == 0:
        return -1e9
    return float((r.log() / arb(10).log()).mid())


def parse(s):
    z = complex(s)
    return acb(z.real, z.imag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=40)
    ap.add_argument('--J', type=int, default=6)
    ap.add_argument('--R', type=int, default=0, help='hoechste x-Ableitung, die gespeichert wird')
    ap.add_argument('--K', type=int, default=6, help='ganzzahlige Hoehen 0..K mit Ableitungen bis J-1')
    ap.add_argument('--prec', type=int, default=800)
    ap.add_argument('--name', default=None)
    ap.add_argument('--budget', type=float, default=1080.0)
    a = ap.parse_args()
    E.setprec(a.prec)
    t0 = time.time()
    pkt = {p: parse(p) for p in PUNKTE}
    werte = {p: [] for p in PUNKTE}          # werte[p][n][r][j]
    ganz = {str(k): [] for k in range(a.K + 1)}  # ganz[k][n][r][j], r<=J-1
    maxrad = -1e9
    g = [E.d_base()] + [{} for _ in range(a.J - 1)]
    grade = []
    for n in range(a.N + 1):
        if time.time() - t0 > a.budget:
            print('Zeitbudget erreicht bei n =', n); break
        # Ableitungen
        ders = [g]
        for r in range(max(a.R, a.J - 1)):
            ders.append([E.deriv(x) if x else {} for x in ders[-1]])
        for p, z in pkt.items():
            row = []
            for r in range(a.R + 1):
                vals = []
                for x in ders[r]:
                    v = E.evaluate(x, z) if x else acb(0)
                    maxrad = max(maxrad, lograd(v))
                    vals.append([v.real.mid().str(60, radius=False), v.imag.mid().str(60, radius=False)])
                row.append(vals)
            werte[p].append(row)
        for k in range(a.K + 1):
            row = []
            for r in range(a.J):
                vals = []
                for x in ders[r]:
                    v = E.evaluate(x, acb(k)) if x else acb(0)
                    maxrad = max(maxrad, lograd(v))
                    vals.append(v.real.mid().str(60, radius=False))
                row.append(vals)
            ganz[str(k)].append(row)
        grade.append([E.degree(x) for x in g])
        print(n, round(time.time() - t0, 1), 's  Grad', grade[-1][-1], ' log10 maxrad %.0f' % maxrad, flush=True)
        if n < a.N:
            g = E.rank_step(g, a.J)
    out = {'N': len(grade) - 1, 'J': a.J, 'R': a.R, 'K': a.K, 'prec': a.prec, 'sekunden': time.time() - t0,
           'log10_max_kugelradius': maxrad, 'grade': grade, 'punkte': PUNKTE, 'werte': werte, 'ganz': ganz}
    name = a.name or f'N{a.N}J{a.J}'
    (HIER / 'data' / f'raenge_{name}.json').write_text(json.dumps(out), encoding='utf-8')
    print('gespeichert', name, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
