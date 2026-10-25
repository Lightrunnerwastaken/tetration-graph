"""D2: T, T', D = dT/db, D' = dT'/db auf einem x-Gitter, Basis als GP-Ausdruck (exakt, z.B. 6/5, 10^4, exp(1)).
Zentrale Differenzen: Basis b +- 10^-HK, Hoehe x +- 10^-EK.  Elastizitaeten: E(x) = b ln b D/T = dlnT/dlnlnb, eps0 = b ln b D'(0)/T'(0).
Aufruf: python d2_elastizitaet.py DPS HK EK GITTER b1 b2 ...   GITTER: 'fein' (-1.95..2, 0.05) oder 'rand' (nur -1, 0)
"""
import sys, json, time, mpmath as mp
from motor import werte, schliessen, HIER

dps = int(sys.argv[1]); HK = int(sys.argv[2]); EK = int(sys.argv[3]); gitter = sys.argv[4]
mp.mp.dps = dps + 10
H = mp.mpf(10) ** (-HK); E = mp.mpf(10) ** (-EK)
if gitter == "fein":
    XS = [mp.mpf(k) / 20 for k in range(-39, 41)]
else:
    XS = [mp.mpf(-1), mp.mpf(0), mp.mpf("-0.5"), mp.mpf("0.5")]


def bwert(bs):
    return mp.mpf(eval(bs.replace("^", "**").replace("exp(", "mp.exp(").replace("/", "*mp.mpf(1)/")))


def xs_gp(xs):
    return [f"sexp({mp.nstr(x, dps + 8)})" for x in xs]


out = {}
for bs in sys.argv[5:]:
    t = time.time(); b = bwert(bs); lb = mp.log(b)
    try:
        Xall = XS + [x + E for x in XS] + [x - E for x in XS]
        V = {}
        for s, tag in ((1, "p"), (-1, "m"), (0, "0")):
            g = f"(({bs})+({s})/10^{HK})" if s else f"({bs})"
            V[tag] = [mp.re(v) for v in werte(g, xs_gp(Xall), dps)]
            schliessen()
        n = len(XS); res = []
        for i, x in enumerate(XS):
            T = V["0"][i]
            T1 = {k: (V[k][n + i] - V[k][2 * n + i]) / (2 * E) for k in V}
            D = (V["p"][i] - V["m"][i]) / (2 * H)
            D1 = (T1["p"] - T1["m"]) / (2 * H)
            res.append({"x": mp.nstr(x, 8), "T": mp.nstr(T, 25), "T1": mp.nstr(T1["0"], 25), "D": mp.nstr(D, 25), "D1": mp.nstr(D1, 20)})
        out[bs] = res
        r0 = [r for r in res if mp.mpf(r["x"]) == 0][0]
        eps0 = b * lb * mp.mpf(r0["D1"]) / mp.mpf(r0["T1"])
        print(f"b={bs}: T'(0)={mp.nstr(mp.mpf(r0['T1']),18)}  D'(0)={mp.nstr(mp.mpf(r0['D1']),15)}  eps0={mp.nstr(eps0,15)}  ({time.time()-t:.0f}s)", flush=True)
    except Exception as e:
        print(bs, "FEHLER", str(e)[-200:], flush=True)
        schliessen()
name = f"d2_{gitter}_{dps}_{HK}_{EK}_" + "_".join(s.replace("/", "o").replace("^", "e").replace("(", "").replace(")", "") for s in sys.argv[5:])
p = HIER / "cache" / f"{name[:150]}.json"
p.write_text(json.dumps(out, indent=0))
print("->", p)
