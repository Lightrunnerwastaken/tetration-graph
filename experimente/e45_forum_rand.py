"""exp.forum_rand — kleinere Forum-Themen nachgerechnet: Ackermann und schnell wachsende Hierarchie, exponentielle Fakultät,
 Grenzfälle bei x → 0, die quadratische Näherung, der Nachbar-Turm x^((1/x)^((1/x)^…)).
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 60
lauf = Lauf("exp.forum_rand", __file__)


# 1. Ackermann–Péter: A(4, n) = 2↑↑(n+3) − 3
def A(m, n):
    if m == 0:
        return n + 1
    if m == 1:
        return n + 2
    if m == 2:
        return 2 * n + 3
    if m == 3:
        return 2 ** (n + 3) - 3
    if n == 0:
        return A(m - 1, 1)
    return A(m - 1, A(m, n - 1))


def turm2(n):
    t = 1
    for _ in range(n):
        t = 2 ** t
    return t


def A_rekursiv(m, n, tiefe=[0]):
    # direkte Rekursion für kleine Fälle (Kontrolle der geschlossenen Formen)
    if m == 0:
        return n + 1
    if n == 0:
        return A_rekursiv(m - 1, 1)
    return A_rekursiv(m - 1, A_rekursiv(m, n - 1))


import sys
sys.setrecursionlimit(100000)
ok_klein = all(A_rekursiv(m, n) == A(m, n) for m in range(4) for n in range(5))
ok_4 = A(4, 0) == turm2(3) - 3 and A(4, 1) == turm2(4) - 3 and A(4, 2) == turm2(5) - 3
lauf.pruefe("A(m, n) für m ≤ 3 geschlossen = Rekursion; A(4, n) = 2↑↑(n+3) − 3 für n = 0, 1, 2", ok_klein and ok_4,
            {"A(4,0)": A(4, 0), "A(4,1)": A(4, 1), "A(4,2) Bits": A(4, 2).bit_length()})


# schnell wachsende Hierarchie: f_0(n) = n+1, f_{k+1}(n) = f_k^[n](n); f_3(n) ≥ 2↑↑n (über log₂-Ketten verglichen)
def f2(n):
    return n * 2 ** n


f3_1 = f2(1)
f3_2 = f2(f2(2))
x = f2(f2(3))                          # f_3(3) = f2(x) = x·2^x, log2 f_3(3) = log2 x + x
log2_f3_3 = mp.log(x, 2) + x
lauf.pruefe("f₃(n) ≥ 2↑↑n für n = 1, 2, 3 (f₃(3) = x·2^x mit x = 402653184·2^402653184)",
            f3_1 >= turm2(1) and f3_2 >= turm2(2) and log2_f3_3 >= turm2(2), {"f3(1)": f3_1, "f3(2)": f3_2, "log2 log2 f3(3)": mp.nstr(mp.log(log2_f3_3, 2), 8)})

# 2. Exponentielle Fakultät: EF(n) = n^EF(n−1); Σ 1/EF(n) ist eine Liouville-Zahl
EF = [1, 1, 2, 9, 262144]              # EF(0) = 1 gesetzt, EF(1..4); EF(5) = 5^262144
s4 = mp.fsum(mp.mpf(1) / e for e in EF[1:])
rest = mp.power(5, -262144)            # 1/EF(5), weitere Glieder vernachlässigbar
summe = s4 + rest
# Liouville: |summe − p/q| < q^{−k} mit q = 2·9·262144 (Nenner der Teilsumme) und grossem k
q = 2 * 9 * 262144
k_liouville = -mp.log(rest) / mp.log(q)
lauf.pruefe("Σ 1/EF(n) = 1.61111492580837673611…, Teilsumme bis n = 4 trifft auf 10^−183000 (Liouville: Exponent k > 10⁴)",
            mp.nstr(summe, 21).startswith("1.6111149258083767361") and k_liouville > 10 ** 4, {"Summe": mp.nstr(summe, 25), "k": mp.nstr(k_liouville, 6)})

# 3. Grenzfälle x → 0+: x↑↑n → 1 (n gerade), → 0 (n ungerade)
def turm(x, n):
    t = mp.mpf(1)
    for _ in range(n):
        t = mp.power(x, t)
    return t


tab = {n: [mp.nstr(turm(mp.mpf(10) ** -e, n), 6) for e in (3, 10, 40)] for n in range(1, 7)}
ok = all(abs(turm(mp.mpf(10) ** -40, n) - (1 if n % 2 == 0 else 0)) < mp.mpf("0.2") for n in range(1, 7))
lauf.pruefe("x↑↑n bei x = 10⁻⁴⁰: gerade n nahe 1, ungerade nahe 0 (langsam: x^x = e^{x ln x})", ok, tab)

# 4. Quadratische Näherung auf (−1, 0]: T(x) = 1 + 2 ln a/(1+ln a)·x − (1 − ln a)/(1+ln a)·x²
def quad(a, x):
    l = mp.log(a)
    return 1 + 2 * l / (1 + l) * x - (1 - l) / (1 + l) * x ** 2


def quad_T(a, x):
    k = int(mp.floor(x))
    y = quad(a, x - k - 1)            # x − k − 1 ∈ [−1, 0)
    for _ in range(k + 1):
        y = mp.power(a, y)
    return y


a = mp.mpf(2)
links = mp.diff(lambda t: quad(a, t), 0)
rechts = mp.diff(lambda t: mp.power(a, quad(a, t - 1)), mp.mpf(0))
zweite_l = mp.diff(lambda t: quad(a, t), 0, 2)
zweite_r = mp.diff(lambda t: mp.power(a, quad(a, t - 1)), mp.mpf(0), 2)
halb = mp.power(2, quad(a, mp.mpf("-0.5")))
lauf.pruefe("quadratische Näherung: T(−1) = 0, C¹ an 0, aber nicht C²; ½-Wert zur Basis 2 = 1.45933 (Kneser 1.45878)",
            abs(quad(a, -1)) < 1e-50 and abs(links - rechts) < 1e-40 and abs(zweite_l - zweite_r) > 0.1 and abs(halb - mp.mpf("1.45933")) < 1e-5,
            {"T(½)": mp.nstr(halb, 10), "T''(0−)": mp.nstr(zweite_l, 6), "T''(0+)": mp.nstr(zweite_r, 6)})

# 5. Nachbar-Turm x^((1/x)^((1/x)^…)) = ln x / W(ln x), konvergent für e^(−1/e) ≤ x ≤ e^e
def nachbar(x, n=4000):
    t = mp.mpf(1)
    for _ in range(n):
        t = mp.power(1 / x, t)
    return mp.power(x, t)


werte = {}
fehler = 0
for x in (mp.mpf("0.8"), mp.mpf(2), mp.e, mp.mpf(10)):
    v = nachbar(x)
    f = mp.log(x) / mp.lambertw(mp.log(x))
    fehler = max(fehler, abs(v - f))
    werte[mp.nstr(x, 4)] = mp.nstr(v, 12)
lauf.pruefe("Nachbar-Turm = ln x / W(ln x) (x = 0.8, 2, e, 10); bei x = e: 1/Ω = 1.7632228343…", fehler < mp.mpf(10) ** -20, werte)
ausserhalb = nachbar(mp.mpf(20), 4001) - nachbar(mp.mpf(20), 4000)
lauf.pruefe("bei x = 20 > e^e konvergiert er nicht (pendelt)", abs(ausserhalb) > 0.1, mp.nstr(ausserhalb, 5))

lauf.daten = {"ef_summe": mp.nstr(summe, 30), "quad_halb_2": mp.nstr(halb, 12), "nachbar": werte}
fig, ax = plt.subplots(1, 1, figsize=(6, 3.6))
xs = [mp.mpf(10) ** (-mp.mpf(i) / 4) for i in range(1, 160)]
for n in range(1, 6):
    ax.semilogx([float(x) for x in xs], [float(turm(x, n)) for x in xs], label=f"x↑↑{n}")
ax.set_xlabel("x"); ax.legend(fontsize=8); ax.set_title("Türme bei x → 0: gerade Höhen → 1, ungerade → 0")
lauf.bild(fig)
lauf.zusammenfassung = "Ackermann A(4,n) = 2↑↑(n+3) − 3, f₃ ≥ 2↑↑n, Σ1/EF(n) = 1.6111149258…, Grenzfälle bei 0, quadratische Näherung C¹ nicht C², Nachbar-Turm = ln x/W(ln x)."
lauf.schreibe()
