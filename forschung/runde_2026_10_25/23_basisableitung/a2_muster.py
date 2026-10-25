"""A2: Vorzeichenmuster, Elastizitaet E = b ln b D/T auf (-1,0), Kriterium max(-1, ln T) < E < 0, und (b > eta)
Vorwaertssumme u0(x) = sum_{k>=0} [r(k) - r(x+k)], r(x) = 1/(b ln^2 b T'(x-1)); Wackeln w = u - u0 (soll 1-periodisch sein).
Aufruf: python a2_muster.py d2_fein_....json"""
import sys, json, mpmath as mp
mp.mp.dps = 40
d = json.load(open(sys.argv[1]))
ETA = mp.exp(mp.exp(-1))


def bwert(bs):
    return mp.mpf(eval(bs.replace("^", "**").replace("exp(", "mp.exp(").replace("/", "*mp.mpf(1)/")))


def inv_tp_bahn(T, T1, b, lb, nmax=5000):
    """Summe_{j>=0} 1/T'(y+j) ab T(y), T'(y) per T(y+1)=b^T(y), T'(y+1)=ln b T(y+1) T'(y)."""
    s = mp.mpf(0); t, tp = T, T1
    for j in range(nmax):
        s += 1 / tp
        if t * lb > 2000:   # naechstes T' riesig
            break
        tn = mp.exp(t * lb); tp = lb * tn * tp; t = tn
        if 1 / tp < mp.mpf(10) ** (-35):
            break
    return s


for bs, rows in d.items():
    b = bwert(bs); lb = mp.log(b)
    R = {mp.nstr(mp.mpf(r["x"]), 6): {k: mp.mpf(v) for k, v in r.items()} for r in rows}
    xs = sorted(R, key=lambda s: mp.mpf(s))
    sig = "".join("+" if R[x]["D"] > mp.mpf(10) ** -20 else ("-" if R[x]["D"] < -mp.mpf(10) ** -20 else "0") for x in xs)
    print(f"b={bs} ({'Kneser' if b > ETA else 'regulaer'})  D auf -1.95..2: {sig}")
    Es = [(mp.mpf(x), b * lb * R[x]["D"] / R[x]["T"], mp.log(R[x]["T"])) for x in xs if -1 < mp.mpf(x) < 0]
    mono = all(e2[1] > e1[1] for e1, e2 in zip(Es, Es[1:]))
    mA = min(e + 1 for _, e, _ in Es); mB = min(e - l for _, e, l in Es)
    # normierte Marge fuer B: (E - lnT)/|lnT|
    mBn = min((e - l) / (-l) for _, e, l in Es)
    print(f"   E auf (-1,0): [{mp.nstr(Es[0][1],8)} .. {mp.nstr(Es[-1][1],8)}], monoton steigend: {mono}; min(E+1)={mp.nstr(mA,6)}; min((E-lnT)/|lnT|)={mp.nstr(mBn,6)}")
    print("   E(-0.75,-0.5,-0.25) =", [mp.nstr(e, 10) for x, e, _ in Es if x in (mp.mpf("-0.75"), mp.mpf("-0.5"), mp.mpf("-0.25"))])
    if b > ETA:
        c = 1 / (b * lb ** 2)
        S0 = c * inv_tp_bahn(R["-1.0"]["T"] if "-1.0" in R else R["-1"]["T"], R["-1.0"]["T1"] if "-1.0" in R else R["-1"]["T1"], b, lb)
        W = {}
        for x in xs:
            xv = mp.mpf(x)
            if xv < mp.mpf("-0.95") or xv > 1:
                continue
            key = mp.nstr(xv - 1, 6)
            if key not in R:
                continue
            Sx = c * inv_tp_bahn(R[key]["T"], R[key]["T1"], b, lb)
            u = R[x]["D"] / R[x]["T1"]
            W[x] = (u, S0 - Sx, u - (S0 - Sx))
        per = max(abs(W[mp.nstr(mp.mpf(x) + 1, 6)][2] - W[x][2]) for x in W if mp.nstr(mp.mpf(x) + 1, 6) in W)
        ws = [W[x][2] for x in W if -1 < mp.mpf(x) <= 0]
        print(f"   Vorwaertssumme: S(0)=sum r(k)={mp.nstr(S0,14)};  Wackeln w=u-u0: Spanne {mp.nstr(max(ws)-min(ws),6)}, max|w|={mp.nstr(max(abs(v) for v in ws),6)}; Periodizitaetsfehler {mp.nstr(per,3)}")
        print("   x, u, u0, w:", [(x, mp.nstr(W[x][0], 8), mp.nstr(W[x][1], 8), mp.nstr(W[x][2], 4)) for x in ("-0.75", "-0.5", "-0.25", "0.5") if x in W])
