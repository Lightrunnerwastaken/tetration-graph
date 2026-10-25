"""D1: Basisableitung D(x) = dT_b(x)/db der Kneser- (b>eta) bzw. regulaeren (b<eta) Tetration auf einem x-Gitter.
Zentrale Differenz in b (Schritt H) und in x (Schritt E), dps DPS.  Ausgabe: D, T, T', u=D/T', rho=b ln b D/T auf (-1,0).
Aufruf: python d1_vorzeichen.py DPS H b1 b2 ...   (b als Dezimalzahl, H z.B. 1e-12)
"""
import sys, json, time, mpmath as mp
from motor import werte, schliessen, HIER

dps = int(sys.argv[1]); H = mp.mpf(sys.argv[2]); mp.mp.dps = dps + 10
E = mp.mpf("1e-12")
XS = [mp.mpf(k) / 20 for k in range(-39, 41)]   # -1.95 .. 2.00
XS = [x for x in XS if x not in (-1, 0)] + [mp.mpf(-1), mp.mpf(0)]
out = {}
for bs in sys.argv[3:]:
    b = mp.mpf(bs); t = time.time()
    try:
        def ev(bb, xs):
            return [mp.re(v) for v in werte(f"({mp.nstr(bb, dps + 5)})", [f"sexp({mp.nstr(x, dps + 5)})" for x in xs], dps)]
        T0 = ev(b, XS)
        Tp = ev(b + H, XS); Tm = ev(b - H, XS); schliessen()
        Xe = [x + E for x in XS] + [x - E for x in XS]
        Te = ev(b, Xe); schliessen()
        n = len(XS); res = []
        for i, x in enumerate(XS):
            D = (Tp[i] - Tm[i]) / (2 * H)
            T1 = (Te[i] - Te[n + i]) / (2 * E)
            res.append([mp.nstr(x, 6), mp.nstr(T0[i], 20), mp.nstr(T1, 20), mp.nstr(D, 20)])
        out[bs] = res
        print(f"b={bs} fertig {time.time() - t:.0f}s", flush=True)
    except Exception as e:
        print(bs, "FEHLER", str(e)[-200:], flush=True)
        schliessen()
p = HIER / "cache" / f"d1_{'_'.join(sys.argv[3:])}.json"
p.write_text(json.dumps(out, indent=0))
print("->", p)
