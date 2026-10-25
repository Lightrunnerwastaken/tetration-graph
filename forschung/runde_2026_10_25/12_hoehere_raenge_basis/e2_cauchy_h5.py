"""Ist die Pentation H_5(b,x) bei festem x analytisch in der Basis b an eta_4 (und eta_5)?
Cauchy-Koeffizienten auf Kreisen um b0 mit zwei Radien; obere Halbebene mit dem Motor, untere per Spiegelung
H_5(conj b) = conj H_5(b) (der Motor scheitert bei Im b < 0 in sexpinit). Stimmen die Koeffizienten beider
Radien ueberein, ist H_5 auf dem groesseren Kreis die Randfunktion EINER analytischen Funktion.
Aufruf: python e2_cauchy_h5.py <b0> <M> <dps> <radius1> <radius2> ..."""
import sys, json, time
import mpmath as mp
import motor

b0 = mp.mpf(sys.argv[1]); M = int(sys.argv[2]); dps = int(sys.argv[3]); rads = [mp.mpf(r) for r in sys.argv[4:]]
mp.mp.dps = dps + 10
XS = ["0.5", "1.5", "-0.5"]
erg = {}
t0 = time.time()
for rho in rads:
    vals = {x: [None] * M for x in XS}
    for j in range(M // 2 + 1):
        b = b0 + rho * mp.expjpi(2 * mp.mpf(j) / M)
        if j == 0 or 2 * j == M:
            b = mp.re(b)
        bs = motor.gpzahl(b, dps + 8)
        w = motor.werte(bs, ["hyper(5,%s)" % x for x in XS], dps)
        for x, v in zip(XS, w):
            vals[x][j] = v
            if 0 < j < M - j:
                vals[x][M - j] = mp.conj(v)
        print(f"rho={rho} j={j} b={mp.nstr(b,8)} H5(0.5)={mp.nstr(w[0],15)} {time.time()-t0:.0f}s", flush=True)
        motor.schliessen()
    for x in XS:
        c = [mp.fsum(vals[x][j] * mp.expjpi(-2 * mp.mpf(j * k) / M) for j in range(M)) / M / rho ** k for k in range(M // 2)]
        erg.setdefault(x, {})[str(rho)] = [mp.nstr(mp.re(ck), 20) for ck in c] + ["im:"] + [mp.nstr(mp.im(ck), 3) for ck in c]
print()
for x in XS:
    print("x =", x)
    r1, r2 = [str(r) for r in rads[:2]]
    for k in range(M // 2):
        a1 = mp.mpf(erg[x][r1][k]); a2 = mp.mpf(erg[x][r2][k])
        print(f"  c_{k}: {mp.nstr(a1,16):>24} {mp.nstr(a2,16):>24}  diff {mp.nstr(abs(a1-a2),3)}")
out = motor.HIER / "daten" / f"cauchy_h5_{sys.argv[1]}.json"
out.write_text(json.dumps({"b0": sys.argv[1], "M": M, "dps": dps, "coeff": erg}, indent=1), encoding="utf-8")
