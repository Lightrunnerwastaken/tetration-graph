"""Abbildungen: Glieder der eps-Reihe (Rang 1 und Halbrang) und Ranggrowth der Ordnungen."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
from halbschritt import lade, folge
r1 = json.loads((HIER / 'data' / 'rang1_J52.json').read_text())
a = json.loads((HIER / 'data' / 'auswertung_N40J8R7.json').read_text())
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for z in r1['zeilen'][:4]:
    t = [abs(complex(float(x), float(y))) for x, y, _ in z['terme_hp']]
    ax[0].semilogy(range(1, len(t) + 1), t, '.-', label='Rang 1, x=%g' % z['x'][0])
for p in ['0.25', '0.5', '0.75', '1.5']:
    t = [abs(complex(float(x), float(y))) for x, y in a['halbrang_jet'][p]['u'][:7]]
    ax[0].semilogy(range(1, 8), t, 'k--', alpha=.5)
ax[0].set_xlabel('Ordnung j'); ax[0].set_ylabel('|Glied der eps-Reihe| bei eps=1')
ax[0].set_title('Rang 1 (farbig) und Halbrang (gestrichelt)'); ax[0].legend(fontsize=8)
d = lade('N40J12')
for j in range(12):
    f = folge(d, '0.5', j)
    ax[1].semilogy(range(len(f)), [abs(complex(v)) + 1e-300 for v in f], label='j=%d' % (j + 1))
ax[1].set_xlabel('ganzzahliger Rang n'); ax[1].set_ylabel('|c_n^(j)(0.5)|'); ax[1].set_ylim(1e-6, 1e125)
ax[1].set_title('Wachstum der Rangjets bei x=0.5'); ax[1].legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(HIER / 'figures' / 'glieder_und_rangwachstum.png', dpi=130)
print('ok')
