"""Etappe 4: schneller Rang-1-Jet von T(f_eps) (O(J^2) Produkte statt O(J^3)).

e(x+1)-e(x) = eps*[2*2^x*E(x) - x - 2 - e(x)],  E = exp(ln2*e) = 2^e,
E_k = (ln2/k) sum_{i=1..k} i e_i E_{k-i}.  Gleiche Auswahl (Klasse A) wie exppoly.antidiff.
Speichert Werte aller Ordnungen an vielen Hoehen (reell und komplex)."""
import sys, json, time, argparse
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import exppoly as E
from flint import acb, arb_poly, arb

def x2(a):  # Multiplikation mit 2^x: Exponent m -> m+1
    return {m + 1: p for m, p in a.items()}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--J', type=int, default=60); ap.add_argument('--prec', type=int, default=2000)
    ap.add_argument('--name', default=None); ap.add_argument('--budget', type=float, default=1100)
    a = ap.parse_args(); E.setprec(a.prec)
    t0 = time.time()
    e = []; Es = [{0: arb_poly([1])}]
    for j in range(1, a.J + 1):
        G = E.scale(x2(Es[j - 1]), 2)
        if j == 1:
            G = E.add(G, {0: arb_poly([-2, -1])})
        else:
            G = E.add(G, e[j - 2], -1)
        # Delta e_j(x) = G(x)  (rechte Seite schon bei x ausgewertet: e(x+1)-e(x)=eps*G(x))
        ej = E.antidiff(G)
        e.append(ej)
        acc = {}
        for i in range(1, j + 1):
            acc = E.add(acc, E.scale(E.mul(e[i - 1], Es[j - i]), i))
        Es.append(E.scale(acc, E.LN2 / j))
        if j % 10 == 0:
            print(j, round(time.time() - t0, 1), 's', flush=True)
        if time.time() - t0 > a.budget:
            print('Budget erreicht bei j =', j); break
    pkt = [k / 16 for k in range(-8, 33)] + [complex(0.5, y) for y in (0.25, 0.5, 0.75, 1.0)] + [complex(0.25, 0.25), complex(0.75, 0.5)]
    werte = {}
    for z in pkt:
        row = []
        for ej in e:
            v = E.evaluate(ej, acb(z.real, z.imag))
            row.append([v.real.mid().str(40, radius=False), v.imag.mid().str(40, radius=False), float((v.rad().log() / arb(10).log()).mid()) if v.rad() != 0 else -999])
        werte[repr(z)] = row
    out = {'J': len(e), 'prec': a.prec, 'sekunden': time.time() - t0, 'werte': werte}
    (HIER / 'data' / f'e4_rang1_{a.name or "J%d" % len(e)}.json').write_text(json.dumps(out), encoding='utf-8')
    print('fertig', len(e), round(time.time() - t0, 1), 's')

if __name__ == '__main__':
    main()
