"""Strang C, Runde 4, c26b: Gilt mein Befund 3a/§3.5 (p_N konvergiert auch im Halbmond D(a0,R0)∖D̄) für grosse N? (Strang B, 22:22)
Π_N f_K aus gemeinsam/lin_PiN_fK_bis_1000.json an Punkten des Halbmonds; Differenzen aufeinanderfolgender N (kein Motor nötig)."""
import json, os, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 70
HIER = os.path.dirname(os.path.abspath(__file__))
PI = json.load(open(os.path.join(HIER, "..", "gemeinsam", "lin_PiN_fK_bis_1000.json"), encoding="utf-8"))
NS = sorted(int(k) for k in PI if k != "beschreibung")
NS = [n for n in NS if n % 50 == 0 or n in (125, 175, 225, 275, 325, 375, 425, 475)]
pkt = [mp.mpf("1.6"), mp.mpf("1.8"), mp.mpc("1.4", "0.9"), mp.mpf("1.84")]
a0 = mp.mpf(PI["1000"]["b"][0]); L = mp.mpc("0.31813150520476413", "1.3372357014306895"); R0 = abs(L - a0)
print("R0 =", mp.nstr(R0, 6), " a0+R0 =", mp.nstr(a0 + R0, 6), " |x| der Punkte:", [mp.nstr(abs(x), 4) for x in pkt],
      " |x−a0|/R0:", [mp.nstr(abs(x - a0) / R0, 4) for x in pkt])
werte = {}
for N in NS:
    c = [mp.mpf(s) for s in PI[str(N)]["b"]]
    werte[N] = [mp.polyval(c[::-1], x) for x in pkt]
print("N | " + " | ".join(f"Π_N f_K({mp.nstr(x, 3)})" for x in pkt) + " | |Π_N − Π_1000| je Punkt")
for N in NS:
    print(f"{N} | " + " | ".join(mp.nstr(v, 12) for v in werte[N]) + " | "
          + " ".join(mp.nstr(abs(werte[N][i] - werte[1000][i]), 3) for i in range(len(pkt))))
