"""D3: dichte Basis-Abtastung: ist b -> T_b(x) bei festem x monoton (steigend auf (-2,-1) und (0,oo), fallend auf (-1,0))?
Basen: rationale Werte von 101/100 bis 1000 (regulaer unter eta, Kneser darueber).  Aufruf: python d3_dicht.py DPS"""
import sys, json, mpmath as mp
from motor import werte, schliessen, HIER
dps = int(sys.argv[1]); mp.mp.dps = dps
XS = ["-1.9", "-1.5", "-1.1", "-0.9", "-0.5", "-0.1", "0.1", "0.5", "0.9", "1.5"]
BS = [f"{n}/100" for n in (101, 103, 106, 110, 115, 120, 125, 130, 135, 140, 142, 143, 144)] + ["1444/1000", "1445/1000", "1446/1000"] + \
     [f"{n}/100" for n in (145, 147, 150, 155, 160, 170, 180, 190, 200, 225, 250, 275, 300, 350, 400, 500, 700, 1000, 1500, 2000, 3000, 5000, 10000, 30000, 100000)]
tab = {}
for bs in BS:
    try:
        tab[bs] = [mp.re(v) for v in werte(f"({bs})", [f"sexp({x})" for x in XS], dps)]
    except Exception as e:
        print(bs, "FEHLER", str(e)[-120:], flush=True)
    schliessen()
ok = [b for b in BS if b in tab]
for j, x in enumerate(XS):
    seq = [tab[b][j] for b in ok]
    dif = [s2 - s1 for s1, s2 in zip(seq, seq[1:])]
    pos = sum(1 for d in dif if d > 0); neg = sum(1 for d in dif if d < 0)
    print(f"x={x:>5}: steigend {pos}, fallend {neg} von {len(dif)};  T von {mp.nstr(seq[0],8)} bis {mp.nstr(seq[-1],8)}")
(HIER / "cache" / f"d3_{dps}.json").write_text(json.dumps({b: [mp.nstr(v, dps) for v in tab[b]] for b in ok}, indent=0))
