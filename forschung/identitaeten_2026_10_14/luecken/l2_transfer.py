"""Transfer (S8/S9), Zahlenprobe zu BERICHT_LUECKEN_ANALYTISCH.md, Abschnitt 5.

Bewiesen wird dort:  a_n − b_n = O(1/(n² ln n)),  a_n − b2_n = O(1/(n³ ln n)),
b_n = −h'(v_n)/(n ln b) + O(1/(n ln² n)),  v_n = (c − ln n)/ln b.
Hier: die skalierten Reste aus den vorhandenen Daten von Etappe 2 (a_n auf 44 Stellen, b_n per Cauchy-Integral).
Aufruf: python l2_transfer.py   (~1 min, kein Motor)
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E2 = HIER.parent / "etappe2"
sys.path.insert(0, str(E2))
from kern import Fix  # noqa: E402

mp.mp.dps = 40
cx = lambda v: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))
erg = {}
for basis, fa, fm, fk in [("e", "a_e_80_0.003_16384.json", "modell_e_80_0.003_16384.json", "korrektur_e.json"),
                          ("2", "a_2_80_0.003_16384.json", "modell_2_80_0.003_16384.json", "korrektur_2.json")]:
    fx = Fix(mp.e if basis == "e" else mp.mpf(2), K=40)
    lb = fx.lb
    a = [cx(v) for v in json.loads((E2 / "daten" / fa).read_text())["a"]]
    md = json.loads((E2 / "daten" / fm).read_text())
    b = [cx(v) for v in md["b"]]
    c = mp.mpc(complex(md["c"].replace(" ", ""))) if "c" in md else None
    A = mp.mpf(md["A"])
    c = mp.log(1j * A / (2 * mp.pi * lb))
    d1 = mp.mpc(complex(json.loads((E2 / "daten" / fk).read_text())["delta1"].replace(" ", "")))
    print(f"Basis {basis}: A = {mp.nstr(A, 15)}, c = {mp.nstr(c, 15)}, δ1 = {mp.nstr(d1, 12)}")
    print("     n   n²ln n|a−b|   n³ln n|a−b2|   |a/T1 − 1|   ln n·|a/T1−1|   n ln²n |b−T1|")
    zeilen = []
    for n in [4, 16, 64, 256, 1024, 2048, 4000]:
        b2 = b[n] - d1 * ((n + 1) * b[n + 1] - 2 * n * b[n] + (n - 1) * b[n - 1])
        v = (c - mp.log(n)) / lb
        hp = mp.diff(fx.h, v)
        T1 = -hp / (n * lb)
        ln = mp.log(n)
        r = [n, n * n * ln * abs(a[n] - b[n]), n ** 3 * ln * abs(a[n] - b2), abs(a[n] / T1 - 1),
             ln * abs(a[n] / T1 - 1), n * ln ** 2 * abs(b[n] - T1)]
        zeilen.append([r[0]] + [mp.nstr(x, 5) for x in r[1:]])
        print(f"{n:6d}  " + "  ".join(f"{mp.nstr(x, 5):>13}" for x in r[1:]))
    erg[basis] = zeilen
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / "l2_transfer.json").write_text(json.dumps(erg, indent=1), encoding="utf-8")
