"""Sprung der Hexation an eta_5: L(x) = Grenzwert der regulaeren Hexation (b -> eta_5-), R(x) = Grenzwert der Kneser-Hexation
(b -> eta_5+), beide per kubischem Einseiten-Fit (Abstaende 0.0025..0.01). Ist der Sprung eine Hoehenverschiebung
R(x) = L(x + theta(x)) mit 1-periodischem theta? Dann muss theta(x) = (R(x)-L(x))/L'(x) bei x und x+1 gleich sein.
Aufruf: python e5_theta_sprung.py <dps>"""
import sys, json, time
import mpmath as mp
import motor

dps = int(sys.argv[1]); mp.mp.dps = dps + 10
ETA5 = mp.mpf("1.737355925151687573101390")
X0 = [mp.mpf(s) for s in ["-0.5", "0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2.5"]]
h = mp.mpf("1e-6")
XS = []
for x in X0: XS += [x - h, x, x + h]
D = [mp.mpf(s) for s in ["0.0025", "0.005", "0.0075", "0.01"]]
t0 = time.time(); W = {}
for seite in (-1, 1):
    for d in D:
        b = motor.gpzahl(ETA5 + seite * d, dps + 8)
        W[(seite, d)] = [mp.re(v) for v in motor.werte(b, ["hyper(6,%s)" % motor.gpzahl(x, dps + 8) for x in XS], dps)]
        motor.schliessen()
        print(seite, mp.nstr(d, 3), f"{time.time()-t0:.0f}s", flush=True)


def grenz(seite, i):
    ts = [seite * d for d in D]; ys = [W[(seite, d)][i] for d in D]
    A = mp.matrix([[t ** k for k in range(4)] for t in ts])
    return mp.lu_solve(A, mp.matrix(ys))[0]


aus = {}
for j, x in enumerate(X0):
    Lm, L0, Lp = [grenz(-1, 3 * j + k) for k in range(3)]
    R0 = grenz(1, 3 * j + 1)
    dL = (Lp - Lm) / (2 * h)
    th = (R0 - L0) / dL
    aus[mp.nstr(x, 4)] = {"L": mp.nstr(L0, 15), "R": mp.nstr(R0, 15), "R-L": mp.nstr(R0 - L0, 6), "L'": mp.nstr(dL, 10), "theta": mp.nstr(th, 6)}
    print(f"x={mp.nstr(x,4):>5}  L={mp.nstr(L0,13)}  R-L={mp.nstr(R0-L0,5):>10}  L'={mp.nstr(dL,8)}  theta={mp.nstr(th,5)}", flush=True)
(motor.HIER / "daten" / f"theta_sprung_eta5_{dps}.json").write_text(json.dumps(aus, indent=1), encoding="utf-8")
