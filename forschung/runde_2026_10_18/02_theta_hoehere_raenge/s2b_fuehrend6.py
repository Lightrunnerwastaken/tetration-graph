"""Rang 6: fuehrender Term a_n ~ -e^{-2 pi i n x*} h'(v_n)/(n ln lam_p), v_n = (c6 - ln n)/ln lam_p (h' per Differenzenquotient).
Aufruf: python s2b_fuehrend6.py"""
import json
from pathlib import Path
import mpmath as mp
import motor
from lib6 import Rang6
HIER = Path(__file__).resolve().parent
R = Rang6("2", 40, start=mp.mpc(-1.99294843593232, 1.44512849663921), rho=mp.mpf("0.2"), rad=mp.mpf("0.04"))
mp.mp.dps = 50
m = json.loads((HIER / "daten" / "modell6.json").read_text())
c6, llp, xs = mp.mpc(m["c6"].replace(" ", "").replace("j", "j")) if False else mp.mpc(complex(m["c6"].replace(" ", ""))), mp.mpf(m["ln_lam_p"]), mp.mpf(m["x_stern"])
c6 = mp.mpc(*[mp.mpf(t) for t in m["c6"].strip("()").replace(" ", "").replace("j", "").replace("+-", "+-").split("+", 1)]) if "+" in m["c6"].strip("()")[1:] else c6
a = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in json.loads((HIER / "daten" / "theta6_40_0.04_512.json").read_text())["a"]]
eps = mp.mpf("1e-10")
out = []
for n in [4, 16, 64, 128]:
    v = (c6 - mp.log(n)) / llp
    hp, hm = R.h6([v + eps, v - eps])[0]
    T = -(hp - hm) / (2 * eps) / (n * llp)
    q = T / (a[n] * mp.expjpi(2 * n * xs))
    print(f"n={n:4d}: T_n / (a_n e^(2 pi i n x*)) = {mp.nstr(q, 8)}")
    out.append([n, str(q)])
(HIER / "daten" / "fuehrend6.json").write_text(json.dumps(out, indent=1))
