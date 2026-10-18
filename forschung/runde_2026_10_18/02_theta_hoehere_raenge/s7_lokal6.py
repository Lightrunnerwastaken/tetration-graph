"""Rang 6: algebraische Singularitaet von theta6 bei x*.
Herleitung: theta6(z) = 2 - z + h(ipent(V)) mit V = hex(z-1) -> -inf;  slog(V) = -2 + 2^V/s1 + ...,  ipent(V) = psing + ipent'(-2) 2^V/s1 + ...
=> theta6(z) = theta* + K (z-x*)^alpha + O(|z-x*|^(2 alpha)),  alpha = ln b / ln lam_p,
   theta* = 2 - x* + h(psing),  K = h'(psing) ipent'(-2) A6^alpha 2^(-pcent) / s1,  s1 = sexp'(-1) = sexp'(0)/ln b.
=> a_n e^{2 pi i n x*} ~ K (i/2pi)^alpha n^(-1-alpha) / Gamma(-alpha).
Aufruf: python s7_lokal6.py"""
import json
from pathlib import Path
import mpmath as mp
import motor
from lib6 import Rang6
HIER = Path(__file__).resolve().parent
DPS = 40
R = Rang6("2", DPS, start=mp.mpc(-1.99294843593232, 1.44512849663921), rho=mp.mpf("0.2"), rad=mp.mpf("0.04"))
mp.mp.dps = DPS + 10
m = json.loads((HIER / "daten" / "modell6.json").read_text())
xs, llp, pcent, A6 = mp.mpf(m["x_stern"]), mp.mpf(m["ln_lam_p"]), mp.mpf(m["pcent"]), mp.mpf(m["A6"])
motor.roh("2", ["pentx(0.5)", "hexx(0.5)"], DPS)
psing = motor.werte("2", ["psing"], DPS)[0]
d = mp.mpf("1e-12")
ip = motor.werte("2", [f"ipentx({motor.gpzahl(-2 + d, DPS + 10)})", f"ipentx({motor.gpzahl(-2 - d, DPS + 10)})", "ipentx(-2)"], DPS)
ipd = (ip[0] - ip[1]) / (2 * d)
print("psing =", mp.nstr(psing, 25), " ipent(-2)+1 =", mp.nstr(ip[2] + 1, 25), " ipent'(-2) =", mp.nstr(ipd, 15))
A4 = mp.mpf("0.8893649546209763728"); s1 = A4 / mp.log(2)
alpha = mp.log(2) / llp
e = mp.mpf("1e-10")
hv, _ = R.h6([psing, psing + e, psing - e])
h0, hd = hv[0], (hv[1] - hv[2]) / (2 * e)
theta_stern = 2 - xs + h0
K = hd * ipd * A6 ** alpha * mp.mpf(2) ** (-pcent) / s1
print("alpha =", mp.nstr(alpha, 20), " theta* =", mp.nstr(theta_stern, 20), " K =", mp.nstr(K, 15))
a = [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in json.loads((HIER / "daten" / "theta6_40_0.04_512.json").read_text())["a"]]
pref = K * (1j / (2 * mp.pi)) ** alpha / mp.gamma(-alpha)
zeilen = []
print(" n    a_n e^{2pi i n x*} / Vorhersage      lokaler Exponent -log2(|a_2n|/|a_n|)")
for n in [5, 10, 20, 40, 80, 160]:
    q = a[n] * mp.expjpi(2 * n * xs) / (pref * mp.mpf(n) ** (-1 - alpha))
    le = -mp.log(abs(a[2 * n]) / abs(a[n]), 2) if 2 * n < 256 else None
    print(f"{n:4d}   {mp.nstr(q, 8):>32}   {mp.nstr(le, 6) if le else ''}")
    zeilen.append([n, str(q), str(le)])
# theta auf x* + iy aus der Fourier-Reihe gegen theta* + K (iy)^alpha
x0 = xs + 3
for y in ["0.08", "0.04", "0.02", "0.01"]:
    yy = mp.mpf(y); z = x0 + 1j * yy
    th = mp.fsum(a[n] * mp.expjpi(2 * n * z) for n in range(len(a)))
    mod = theta_stern + K * (1j * yy) ** alpha
    print(f"y={y}: theta={mp.nstr(th, 12)}  theta-theta*={mp.nstr(th-theta_stern,6)}  Rest nach K-Term={mp.nstr(th-mod,4)}  Rest/y^(2alpha)={mp.nstr(abs(th-mod)/yy**(2*alpha),4)}")
(HIER / "daten" / "lokal6.json").write_text(json.dumps({"alpha": str(alpha), "theta_stern": str(theta_stern), "K": str(K), "psing": str(psing), "zeilen": zeilen}, indent=1))
