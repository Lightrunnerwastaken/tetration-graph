"""Etappe 6: einheitlich Kneser-artige Raenge H4, H5, H6, ... zur Basis 2 (Cauchy-Verfahren, e6_superfunktion).
Speichert die Mittellinienwerte (h, c, L) jeder Stufe, damit spaetere Laeufe darauf aufbauen.
Aufruf: python e6_raenge.py Rmax A N"""
import sys, json, time, pickle
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import numpy as np
from e6_superfunktion import Superfunktion, Exp2, fixpunkt, LN2
Rmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
Rvon = int(sys.argv[4]) if len(sys.argv) > 4 else 4  # Stufen < Rvon aus Dateien laden
A = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
N = int(sys.argv[3]) if len(sys.argv) > 3 else 500
t0 = time.time()
unten = Exp2(); L = fixpunkt(unten.f, 0.82 + 1.57j)
stufen = {}
jf = HIER / 'data' / ('e6_raenge_A%g_N%d.json' % (A, N))
if Rvon > 4 and jf.exists():
    stufen = {int(k): v for k, v in json.loads(jf.read_text()).items()}
for r in range(4, Rmax + 1):
    pk = HIER / 'data' / ('e6_H%d_A%g_N%d.pkl' % (r, A, N))
    if r < Rvon:
        with open(pk, 'rb') as fh:
            S = Superfunktion.laden(unten, pickle.load(fh), name='H%d' % r)
        unten = S; L = fixpunkt(S, L); print('geladen H%d' % r, flush=True); continue
    S = Superfunktion(unten, L, A=A, N=N, name='H%d' % r)
    h0 = None
    if r > 4:
        Lp = unten.L
        h0 = unten.h + np.where(unten.y >= 0, L - Lp, np.conj(L - Lp)) * np.tanh(np.abs(unten.y) / 2)
    res = S.loese(h0=h0, iters=int(sys.argv[5]) if len(sys.argv) > 5 else 3000, lam=0.3, tol=1e-13)
    # Kontrollen
    x = np.array([0.25, 0.5, 0.75])
    # echte Kontrolle: direkte Cauchy-Auswertung bei w+1 gegen f(F(w)) im Rechteck
    wv = np.array([-0.3, -0.2, -0.3 + 0.3j]) + 0j
    fg = np.max(np.abs(S._F(wv + 1) - unten.f(S._F(wv))))
    Lneu = fixpunkt(S, L)
    lam = S.ableitung(Lneu)[0]
    stufen[r] = {'L_unten': [L.real, L.imag], 'residuum': float(res), 'iter': S.iter, 'c': S.c,
                 'werte_0.5_1.5': [[complex(v).real, complex(v).imag] for v in S(np.array([0.5, 1.5, 0.25, 0.75, -0.5]))],
                 'funktionalgleichung': float(fg), 'fixpunkt': [Lneu.real, Lneu.imag], 'multiplikator': [lam.real, lam.imag]}
    print('H%d: Res %.1e (%d It.), c=%.6f, H(0.5)=%s, H(1.5)=%s, |S(x+1)-f(S(x))|=%.1e, Fixpunkt %s, Mult %s, %.0fs' % (
        r, res, S.iter, S.c, S(0.5)[0], S(1.5)[0], fg, Lneu, lam, time.time() - t0), flush=True)
    with open(HIER / 'data' / ('e6_H%d_A%g_N%d.pkl' % (r, A, N)), 'wb') as fh:
        pickle.dump(S.speichern(), fh)
    (HIER / 'data' / ('e6_raenge_A%g_N%d.json' % (A, N))).write_text(json.dumps(stufen, indent=1), encoding='utf-8')
    unten = S; L = Lneu
