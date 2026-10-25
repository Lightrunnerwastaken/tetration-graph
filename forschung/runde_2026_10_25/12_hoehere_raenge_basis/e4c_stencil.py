"""Wie e4b, aber mit frei gewaehlten Abstaenden (Stencil) -- trennt Abbruchfehler des Fits von einem echten Defekt.
Einseitige Polynomfits an eta_5: stimmen Wert und Ableitungen von links und rechts ueberein?
Aufruf: python e4b_kontrolle.py <dps> <zentrum> <rang> <abstaende,durch,komma>"""
import sys, json, time
import mpmath as mp
import motor

dps = int(sys.argv[1]); mp.mp.dps = dps + 10
ETA5 = mp.mpf(sys.argv[2]); RANG = int(sys.argv[3])
XS = ["0.5", "1.5"]
D = [mp.mpf(s) for s in sys.argv[4].split(",")]
t0 = time.time()
werte = {}
for seite in (-1, 1):
    for d in D:
        b = ETA5 + seite * d
        try:
            w = motor.werte(motor.gpzahl(b, dps + 8), ["hyper(%d,%s)" % (RANG, "%s") % x for x in XS], dps)
            werte[(seite, d)] = [mp.re(v) for v in w]
            print(seite, mp.nstr(d, 3), [mp.nstr(mp.re(v), 18) for v in w], f"{time.time()-t0:.0f}s", flush=True)
        except Exception as e:
            print(seite, mp.nstr(d, 3), "FEHLER", str(e)[:200], flush=True)
        motor.schliessen()
erg = {}
for i, x in enumerate(XS):
    for seite in (-1, 1):
        pts = [(seite * d, werte[(seite, d)][i]) for d in D if (seite, d) in werte]
        if len(pts) < 4: continue
        # kubisches Interpolationspolynom in t = b - eta5, Wert und Ableitungen bei t = 0
        ts = [p[0] for p in pts]; ys = [p[1] for p in pts]
        A = mp.matrix([[t ** k for k in range(4)] for t in ts]); c = mp.lu_solve(A, mp.matrix(ys))
        erg[(x, seite)] = [c[0], c[1], 2 * c[2]]
    if (x, -1) in erg and (x, 1) in erg:
        L, R = erg[(x, -1)], erg[(x, 1)]
        print(f"x={x}: Wert links {mp.nstr(L[0],14)} rechts {mp.nstr(R[0],14)} diff {mp.nstr(L[0]-R[0],3)}; "
              f"Ableitung {mp.nstr(L[1],10)} / {mp.nstr(R[1],10)} diff {mp.nstr(L[1]-R[1],3)}; zweite {mp.nstr(L[2],8)} / {mp.nstr(R[2],8)}")
(motor.HIER / "daten" / f"kontrolle_r{RANG}_{sys.argv[2]}_{sys.argv[4]}_{dps}.json").write_text(json.dumps(
    {f"{s}|{mp.nstr(d,3)}": [mp.nstr(v, dps) for v in w] for (s, d), w in werte.items()}, indent=1), encoding="utf-8")
