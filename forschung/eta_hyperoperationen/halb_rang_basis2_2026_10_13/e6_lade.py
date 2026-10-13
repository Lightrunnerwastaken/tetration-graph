"""Etappe 6: geloeste einheitliche Raenge aus data/e6_H*.pkl laden."""
import pickle
from pathlib import Path
import numpy as np
from e6_superfunktion import Superfunktion, Exp2
HIER = Path(__file__).resolve().parent

def lade_raenge(Rmax, A=30, N=500):
    unten = Exp2(); out = {3: unten}
    for r in range(4, Rmax + 1):
        pk = HIER / 'data' / ('e6_H%d_A%g_N%d.pkl' % (r, A, N))
        if not pk.exists():
            break
        with open(pk, 'rb') as fh:
            S = Superfunktion.laden(unten, pickle.load(fh), name='H%d' % r)
        out[r] = S; unten = S
    return out

def H(raenge, r, z):
    z = np.atleast_1d(np.asarray(z, dtype=complex))
    if r == 3:
        return np.exp(np.log(2.0) * z)
    return raenge[r](z)
