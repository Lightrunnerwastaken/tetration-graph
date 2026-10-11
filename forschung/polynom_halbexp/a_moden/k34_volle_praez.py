"""Strang A: Newton bei grossem N OHNE Zwischenrundung, in Etappen.
  python k34_volle_praez.py N pi        → Π_N f_K in voller Präzision speichern (daten/voll_N.json, alle Stellen)
  python k34_volle_praez.py N schritt   → ein Newton-Schritt vom letzten Iterat (volle Stellen), Korrektur protokolliert
  python k34_volle_praez.py N bericht   → |p − Π| usw."""
import json, sys, time
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from flint import ctx
from lin_kern import DATEN, kneser_koeff, newton_N, pi_N
N, modus = int(sys.argv[1]), sys.argv[2]
prec = int(0.72 * N) + 200
stellen = int(prec * 0.30103) + 30
mp.mp.dps = 60
L = -mp.lambertw(-1); RL = abs(L)
pfad = DATEN / f"voll_{N}.json"
t0 = time.time()
if modus == "pi":
    b = pi_N(kneser_koeff(N + 1), N, prec)
    ctx.prec = prec
    json.dump({"PiN": [x.mid().str(stellen, radius=False) for x in b], "iterate": [], "korr": []}, open(pfad, "w"))
    print(f"N = {N}: Π_N f_K gespeichert ({time.time() - t0:.0f} s)")
elif modus == "schritt":
    d = json.load(open(pfad))
    start = d["iterate"][-1] if d["iterate"] else d["PiN"]
    p, korr = newton_N(start, N, prec, schritte=1, stellen=stellen)
    d["iterate"].append(p); d["korr"].append(mp.nstr(korr[0], 4))
    json.dump(d, open(pfad, "w"))
    print(f"N = {N}: voller Newton-Schritt {len(d['korr'])}: Korrektur {d['korr'][-1]} ({time.time() - t0:.0f} s)")
else:
    d = json.load(open(pfad))
    b = [mp.mpf(x) for x in d["PiN"]]; p = [mp.mpf(x) for x in d["iterate"][-1]]
    print(f"N = {N}: Korrekturen {d['korr']}; max_k |p − Π|·|L|^k = {mp.nstr(max(abs(x - y) * RL ** k for k, (x, y) in enumerate(zip(p, b))), 4)}; "
          f"|p − Π|(L) = {mp.nstr(abs(mp.fsum((x - y) * L ** k for k, (x, y) in enumerate(zip(p, b)))), 4)}")
