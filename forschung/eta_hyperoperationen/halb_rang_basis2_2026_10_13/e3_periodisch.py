"""Etappe 3: Form der wachsenden Moden. Fit A(x)=a*sin(2 pi x)+b*(1-cos(2 pi x)) an x=0.25, 0.5;
Vorhersage an 0.75, -0.5, 0.5+0.5i, 0.25+0.25i, 0.5+1i. (a) Rangmode j=10, n=40 (N40J12);
(b) spaete Glieder j=36 bei Rang 1 und 2 (e3_raenge_J36)."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from halbschritt import lade, folge
mp.mp.dps = 50
def fit(val):
    r = val['0.5']
    A = {p: v / r for p, v in val.items()}
    b = mp.mpf(1) / 2  # A(0.5)=2b=1
    a = mp.re(A['0.25']) - b
    pred = lambda z: a * mp.sin(2 * mp.pi * z) + b * (1 - mp.cos(2 * mp.pi * z))
    return a, {p: (A[p], pred(mp.mpc(complex(p)))) for p in A}
res = {}
d = lade('N40J12')
for j in [9, 11]:
    val = {p: folge(d, p, j)[40] for p in ['0.25', '0.5', '0.75', '-0.5', '0.5+0.5j', '0.25+0.25j', '0.5+1j', '1.5', '1.25']}
    a, tab = fit(val)
    print('Rangmode j=%d, n=40: a=%s' % (j + 1, mp.nstr(a, 5)))
    for p, (ist, soll) in tab.items(): print('   %-11s ist %-28s Modell %s' % (p, mp.nstr(ist, 5), mp.nstr(soll, 5)))
    res['rangmode_j%d' % (j + 1)] = {p: [str(ist), str(soll)] for p, (ist, soll) in tab.items()}
e = json.loads((HIER / 'data' / 'e3_raenge_J36.json').read_text())
for n in ['1', '2']:
    val = {p: mp.mpc(mp.mpf(e['raenge'][n][p][35][0]), mp.mpf(e['raenge'][n][p][35][1])) for p in e['punkte'] if p != '3'}
    a, tab = fit(val)
    print('spaete Glieder Rang %s, j=36: a=%s' % (n, mp.nstr(a, 5)))
    for p, (ist, soll) in tab.items(): print('   %-11s ist %-28s Modell %s' % (p, mp.nstr(ist, 5), mp.nstr(soll, 5)))
    res['spaet_rang%s' % n] = {p: [str(ist), str(soll)] for p, (ist, soll) in tab.items()}
(HIER / 'data' / 'e3_periodisch.json').write_text(json.dumps(res, indent=1), encoding='utf-8')
