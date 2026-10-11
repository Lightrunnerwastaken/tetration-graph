"""Monodromie der q-Abel-Koordinate W = α_q (Fork: knwalk, stetig längs einer Punktfolge) um den anziehenden Fixpunkt p.
 Vergleich: ΔW ∈ ℤ + Πℤ (Π = 2πi/ln P'(q)) bedeutet „Kopie“ (dieselbe Lösung), sonst Kandidat für eine neue verschiebungsfeste reelle Kurve.
 Zum Vergleich ein Umlauf um einen Punkt ohne Singularität und einer um q selbst (muss Π geben).
 Fork nur gelesen. Aufruf: umlauf_w.py
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

mp.mp.dps = 40
t0 = time.time()
P_FIX = complex(1.7937303432391165676, 1.5736302435118265996)
Q = complex(-1.992948435932322332933652, 1.445128496639205463533289)
PI_ = 2j * np.pi / np.log(complex(1.368932747339784991147632, 2.488144667623474107366636))
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)


def umlauf(mitte, r, n=240, name=""):
    pts = [mitte + r * np.exp(1j * f) for f in np.linspace(0, 2 * np.pi, n + 1)]
    liste = "[" + ",".join(f"{w.real:.20f}+I*{w.imag:.20f}" for w in pts) + "]"
    v = gp.eval_batch("2", [f"if(1, kuW=iferr(knwalk({liste}),E,0); if(type(kuW)==\"t_VEC\",1,0))"]
                      + [f"if(type(kuW)==\"t_VEC\", kuW[{j + 1}], 10^9)" for j in range(len(pts))])
    if v[0].real != 1:
        print(f"{name}: knwalk gescheitert [{time.time() - t0:.0f}s]", flush=True)
        return None
    W = np.array([complex(x) for x in v[1:]])
    sprung = np.max(np.abs(np.diff(W)))
    dW = W[-1] - W[0]
    # Zerlegung dW = a + bΠ
    b = dW.imag / PI_.imag
    a = (dW - b * PI_).real
    print(f"{name}: Mitte {mitte:.4f}, r = {r}: ΔW = {dW:.8f} = {a:.6f} + {b:.6f}·Π; grösster Schritt {sprung:.3f}; "
          f"s·Im W von {W.imag.min():.3f} bis {W.imag.max():.3f} [{time.time() - t0:.0f}s]", flush=True)
    return W


try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", "kns"])
    print(f"Π = {PI_:.8f}, s = {float(k[1].real):.0f}", flush=True)
    umlauf(Q, 0.2, name="Kontrolle um q")
    umlauf(complex(-1.2, 0.9), 0.2, name="Kontrolle ohne Singularität")
    for r in (0.15, 0.3, 0.5):
        umlauf(P_FIX, r, name="um p")
finally:
    gp.close()
