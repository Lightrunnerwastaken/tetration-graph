"""A1: Auswertung der D1-Daten: Vorzeichenmuster von D, Kriterium rho = b ln b D/T auf (-1,0):  max(-1, ln T) < rho < 0.
Aufruf: python a1_auswertung.py datei.json"""
import sys, json, mpmath as mp
mp.mp.dps = 30
d = json.load(open(sys.argv[1]))
for bs, rows in d.items():
    b = mp.mpf(bs); lb = mp.log(b)
    rows = sorted(rows, key=lambda r: mp.mpf(r[0]))
    sig = "".join("+" if mp.mpf(r[3]) > 1e-25 else ("-" if mp.mpf(r[3]) < -1e-25 else "0") for r in rows)
    # Kriterium auf (-1,0)
    marg_A = mp.inf; marg_B = mp.inf; rmin = mp.inf; rmax = -mp.inf
    for r in rows:
        x = mp.mpf(r[0])
        if -1 < x < 0:
            T = mp.mpf(r[1]); D = mp.mpf(r[3]); rho = b * lb * D / T
            rmin = min(rmin, rho); rmax = max(rmax, rho)
            marg_A = min(marg_A, rho + 1); marg_B = min(marg_B, rho - mp.log(T))
    # Nullstellen von D (Vorzeichenwechsel) ausser -1, 0
    nulls = []
    for r1, r2 in zip(rows, rows[1:]):
        D1, D2 = mp.mpf(r1[3]), mp.mpf(r2[3])
        if D1 * D2 < 0:
            nulls.append(r1[0] + ".." + r2[0])
    print(f"b={bs:>8}  x=-1.95..2: {sig}")
    print(f"     rho in [{mp.nstr(rmin,6)}, {mp.nstr(rmax,6)}]  min(rho+1)={mp.nstr(marg_A,6)}  min(rho-lnT)={mp.nstr(marg_B,6)}  Wechsel: {nulls}")
    for r in rows:
        if r[0] in ("-1.5", "-0.5", "0.5", "1.5"):
            print(f"     x={r[0]:>5} T={mp.nstr(mp.mpf(r[1]),12)} D={mp.nstr(mp.mpf(r[3]),15)} u=D/T'={mp.nstr(mp.mpf(r[3])/mp.mpf(r[2]),12)}")
