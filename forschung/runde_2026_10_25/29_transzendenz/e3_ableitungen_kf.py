# e3: Kettenbrueche der Ableitungs-Konstanten (sexp'(0), slog'(0), ...) auf 150 Stellen.
# Daten: cache/konstanten_200.json (Kopie aus forschung/identitaeten_2026_10_14, 160 sichere Stellen);
# Gegenprobe: dieselben Konstanten aus konstanten_160.json (andere Genauigkeitsstufe) muessen denselben Kettenbruch geben.
import json, math, sys
from fractions import Fraction
sys.path.insert(0, '.')
from e1_kf_statistik import stats
import e1_kf_statistik as E

E.D = 150
d200 = json.load(open('cache/konstanten_200.json'))
d160 = json.load(open('D:/Users/janis/Documents/Matura/tetration-graph/forschung/identitaeten_2026_10_14/daten/konstanten_160.json'))
keys = ['e.taylor0_a1', 'e.taylor0_a2', 'e.slog\'(0)', 'e.slog\'\'(0)', 'e.T\'\'(0)', 'e.int_0^1_GL140', 'e.int_-1^0_GL140',
        'e.expfak_0!2=F\'(1)_k6', '2.taylor0_a1', '2.int_0^1_GL140']
for k in keys:
    if k not in d200:
        print('fehlt', k); continue
    a = E.cf_interval(str(d200[k]))
    b = E.cf_interval(str(d160.get(k, d200[k]))) if k in d160 else None
    same = (b is not None) and (a[:min(len(a), len(b)) - 2] == b[:min(len(a), len(b)) - 2])
    st = stats(a)
    print('%-26s n=%3d  khin=%.3f  levy=%.4f  mu_max=%.4f  a_max=%d  Gegenprobe(160)=%s  Anfang=%s' % (
        k, st['n'], st['khin'], st['levy'], st['mu_max'], st['a_max'], same, a[:12]))
