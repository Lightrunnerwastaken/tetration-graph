"""A3: Vorwaertssumme S_b(0) = sum_{k>=0} r(k) = 1/(b ln^2 b) * sum_{j>=-1} 1/T'(j)  allein aus T'(0) und den ganzen Tuermen;
Vergleich mit dem Kanalgesetz S ~ (pi/2) sqrt(2 eta/e) (b-eta)^(-3/2).   Aufruf: python a3_kanal.py d2_rand_....json [weitere]"""
import sys, json, mpmath as mp
mp.mp.dps = 40
ETA = mp.exp(mp.exp(-1)); K = mp.pi / 2 * mp.sqrt(2 * ETA / mp.e)


def bwert(bs):
    return mp.mpf(eval(bs.replace("^", "**").replace("exp(", "mp.exp(").replace("/", "*mp.mpf(1)/")))


print(f"Kanalkonstante (pi/2)sqrt(2 eta/e) = {mp.nstr(K, 15)}")
for f in sys.argv[1:]:
    d = json.load(open(f))
    for bs, rows in d.items():
        b = bwert(bs); lb = mp.log(b)
        if b <= ETA:
            continue
        r0 = [r for r in rows if mp.mpf(r["x"]) == 0][0]
        tp0 = mp.mpf(r0["T1"])
        s = lb / tp0          # 1/T'(-1) = ln b / T'(0)
        t, tp = mp.mpf(1), tp0
        for j in range(10 ** 6):
            s += 1 / tp
            if t * lb > 3000:
                break
            t = mp.exp(t * lb); tp = lb * t * tp
        S = s / (b * lb ** 2); dl = b - ETA
        print(f"b={bs}: delta={mp.nstr(dl,6)}  S(0)={mp.nstr(S,14)}  S*delta^1.5={mp.nstr(S*dl**1.5,10)}  Verhaeltnis zu K: {mp.nstr(S*dl**1.5/K,10)}  (Schritte {j})")
