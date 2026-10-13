"""Etappe 4: Zerlegung der spaeten Glieder in theta-Moden und Lage der Borel-Singularitaeten.

Ansatz fuer grosse j:  c_j(x) ~ K'(x) * [g_j e^{2 pi i x} + conj(g_j) e^{-2 pi i x} - 2 Re g_j]
(1-periodisch, null an ganzen Hoehen; K' = Kneser-Ableitung als Naeherung fuer S').
g_j aus x=1/4 und x=3/4 (Projektion), Kontrolle an x=1/2, 1/8, 5/8.
Quotienten g_{j+1}/g_j ~ (j+beta)/A  ->  A aus der Differenz aufeinanderfolgender Quotienten."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e4_theta import lade_reihe, lade_kneser
mp.mp.dps = 40
J, W = lade_reihe(sys.argv[1] if len(sys.argv) > 1 else 'e4_rang1_J84.json')
K = lade_kneser()
h = mp.mpf(1) / 16
Kd = lambda x: (-K[x + 2 * h] + 8 * K[x + h] - 8 * K[x - h] + K[x - 2 * h]) / (12 * h)
def P(x, j):
    return mp.re(W[x][j]) / mp.re(Kd(x))
out = {'g': [], 'kontrolle': []}
for j in range(J):
    u = P(0.25, j); v = P(0.75, j)        # u = a + b, v = -a + b  fuer a sin + b(1-cos)
    a = (u - v) / 2; b = (u + v) / 2
    pred = lambda x: a * mp.sin(2 * mp.pi * x) + b * (1 - mp.cos(2 * mp.pi * x))
    k = [abs(P(x, j) - pred(x)) / (abs(a) + abs(b)) for x in (0.5, 0.125, 0.625, 0.375)]
    g = mp.mpc(-b / 2, -a / 2)  # a sin + b(1-cos) = b - 2 Re(g e^{2pi i x}) mit g = (b + i a)/2 ... Vorzeichen egal fuer |A|
    out['g'].append([str(a), str(b)]); out['kontrolle'].append([float(z) for z in k])
    if j % 4 == 3 or j > J - 6:
        print('j=%2d a=%-12s b=%-12s rel. Modellfehler (0.5,0.125,0.625,0.375): %s' % (j + 1, mp.nstr(a, 5), mp.nstr(b, 5), ' '.join('%.1e' % float(z) for z in k)))
# Borel-Lage aus g_j = (b + i a)/2: Quotienten
gs = [mp.mpc(mp.mpf(b), mp.mpf(a)) / 2 for a, b in out['g']]
print('Quotienten q_j = g_j/g_{j-1} und A = 1/(q_j - q_{j-1}):')
As = []
for j in range(40, J):
    q = gs[j] / gs[j - 1]; q1 = gs[j - 1] / gs[j - 2]
    A = 1 / (q - q1); As.append(A)
    print('  j=%d q=%s  A=%s  |A|=%s' % (j + 1, mp.nstr(q, 5), mp.nstr(A, 5), mp.nstr(abs(A), 5)))
out['A'] = [[float(mp.re(z)), float(mp.im(z))] for z in As]
(HIER / 'data' / 'e4_stokes_folge.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
