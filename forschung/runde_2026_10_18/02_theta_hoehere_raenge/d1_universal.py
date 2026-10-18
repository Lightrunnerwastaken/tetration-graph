"""Frage 3: Ist |a_n(e)| ~ |a_n(2)| ein Gesetz?  Fuehrender Term (kn.satz.theta_koeffizienten):
   a_n ~ T_b(n) = -h_b'(v_n) / (n ln b),  v_n = (c_b - ln n)/ln b,  c_b = ln(A_b/(2 pi ln b)) + i pi/2.
Kettenregel h'(v) = h'(log_b v)/(v ln b) fuer grosse |v|.  Vergleich mit den gemessenen a_n (etappe2, 44 Stellen).
Aufruf: python d1_universal.py"""
import json, sys
from pathlib import Path
import mpmath as mp
E2 = Path(__file__).resolve().parents[2] / "identitaeten_2026_10_14" / "etappe2"
sys.path.insert(0, str(E2))
from kern import Fix  # nur Rechenklasse (kein Motor)
mp.mp.dps = 40
A = {"e": mp.mpf("1.0917673512583209918"), "2": mp.mpf("0.8893649546209763728")}
fx = {"e": Fix(mp.e, K=40), "2": Fix(mp.mpf(2), K=40)}

def hprime(b, v):
    f = fx[b]; acc = mp.mpc(1)
    while abs(v) > 30:
        acc /= v * f.lb
        v = mp.log(v) / f.lb
    return acc * mp.diff(f.h, v)

def T(b, lnn):
    f = fx[b]
    c = mp.log(A[b] / (2 * mp.pi * f.lb)) + 1j * mp.pi / 2
    return -hprime(b, (c - lnn) / f.lb) / f.lb      # = n * a_n (fuehrend)

def lade(b):
    d = json.loads((E2 / "daten" / f"a_{b}_80_0.003_16384.json").read_text())
    return [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in d["a"]]

out = {"gemessen": [], "fuehrend": []}
a = {b: lade(b) for b in "e2"}
print(" n     |a_n(e)|      |a_n(2)|     Verh. gemessen  Verh. fuehrend  arg e   arg 2")
for n in [1, 2, 4, 8, 16, 64, 256, 1024, 2048, 4096]:
    re_ = abs(a["e"][n]) / abs(a["2"][n])
    Te, T2 = T("e", mp.log(n)), T("2", mp.log(n))
    rf = abs(Te) / abs(T2)
    print("%5d  %.5e  %.5e   %.6f       %.6f     %.4f %.4f" % (n, abs(a["e"][n]), abs(a["2"][n]), re_, rf,
          mp.arg(a["e"][n]), mp.arg(a["2"][n])))
    out["gemessen"].append([n, float(abs(a["e"][n])), float(abs(a["2"][n])), float(re_), float(rf)])
print("\nfuehrender Term fuer grosse n:  ln ln n | |T_e|/|T_2| | n ln n |a_n(e)| | n ln n |a_n(2)| | arg e | arg 2")
for lnln in ["1", "2", "2.1178", "3", "4", "5", "7", "10", "15", "20", "30", "50", "100", "300", "1000", "1e4", "1e5", "1e6", "3814279.1", "1e7"]:
    t = mp.exp(mp.mpf(lnln))
    Te, T2 = T("e", t), T("2", t)
    r = abs(Te) / abs(T2)
    print("%10s  %.6f  %.5e  %.5e  %+.4f %+.4f" % (lnln, r, abs(Te) * t, abs(T2) * t, mp.arg(Te), mp.arg(T2)))
    out["fuehrend"].append([lnln, float(r), float(abs(Te) * t), float(abs(T2) * t), float(mp.arg(Te)), float(mp.arg(T2))])
Path("daten/d1_universal.json").write_text(json.dumps(out, indent=1))
