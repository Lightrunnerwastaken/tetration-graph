"""Summen bei eps=1: Tetration (Rang 1, 52 Ordnungen) gegen Kneser, H_3.5 und H_4.5 (Ordnungen 1..7).

Abbruch nach der kleinsten Ordnung ('optimaler Abbruch'); als Unsicherheit wird der Betrag
des kleinsten Glieds angegeben (heuristische Fehlerskala einer divergenten asymptotischen Reihe,
kein Beweis).
"""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
mp.mp.dps = 50


def opt_abbruch(x, terme, jmax=None):
    """Partialsumme bis vor das kleinste Glied (ab Ordnung 3 gesucht)."""
    T = terme[:jmax] if jmax else terme
    absw = [abs(t) for t in T]
    jstar = min(range(2, len(T)), key=lambda j: absw[j])
    s = x + 1 + sum(T[:jstar])
    return s, absw[jstar], jstar + 1


def main():
    out = {}
    r1 = json.loads((HIER / 'data' / 'rang1_J52.json').read_text())
    kn = {complex(*z['x']): mp.mpc(*z['kneser']) for z in r1['zeilen']}
    print('Rang 1 (Tetration), eps=1, optimaler Abbruch gegen Kneser:')
    t1 = {}
    for z in r1['zeilen'][:10]:
        x = complex(*z['x'])
        c = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b, _ in z['terme_hp']]
        s, kl, jst = opt_abbruch(mp.mpc(x), c)
        ref = kn[x]
        # Band der Partialsummen Ordnung 4..12
        ps = [mp.mpc(x) + 1 + sum(c[:j]) for j in range(4, 13)]
        lo = min(mp.re(p) for p in ps); hi = max(mp.re(p) for p in ps)
        t1[str(x)] = {'summe': [float(mp.re(s)), float(mp.im(s))], 'kleinstes_glied': float(kl), 'ordnung': jst,
                      'kneser': [float(mp.re(ref)), float(mp.im(ref))], 'abstand': float(abs(s - ref)),
                      'band_re_4_12': [float(lo), float(hi)]}
        print('  x=%-10s S_opt=%s  Kneser=%s  |S-K|=%.1e  kleinstes Glied %.1e (j=%d)  Band Re[4..12]=[%.6f, %.6f]' % (
            x, mp.nstr(s, 9), mp.nstr(ref, 9), abs(s - ref), kl, jst, lo, hi))
    out['rang1'] = t1
    a = json.loads((HIER / 'data' / 'auswertung_N40J8R7.json').read_text())
    for key, name in [('u', 'H_3.5'), ('u32', 'H_4.5')]:
        print(name + ' (Ordnungen 1..7), eps=1:')
        tab = {}
        for p, v in a['halbrang_jet'].items():
            x = mp.mpc(complex(p))
            terme = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re, im in v[key]][:7]
            s, kl, jst = opt_abbruch(x, terme)
            s7 = x + 1 + sum(terme)
            tab[p] = {'summe_opt': [float(mp.re(s)), float(mp.im(s))], 'kleinstes_glied': float(kl), 'ordnung': jst,
                      'summe_1_7': [float(mp.re(s7)), float(mp.im(s7))], 'glieder': [float(abs(t)) for t in terme]}
            print('  x=%-11s S_opt=%-28s kleinstes Glied %.1e (j=%d)  |Glieder|: %s' % (
                p, mp.nstr(s, 8), kl, jst, ' '.join('%.0e' % float(abs(t)) for t in terme)))
        out[name] = tab
    (HIER / 'data' / 'summen.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
