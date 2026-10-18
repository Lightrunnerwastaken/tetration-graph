"""E1: Gegenproben fuer tower_mod und exakte Bewertungsformel
   v_p(a^^H - alpha_p(a)) = V_p(a) (H - t_p(a)) fuer H >= t_p(a)   (p ungerade, p ∤ a)
   v_2(a^^H - alpha_2(a)) = v_2(a-1) + H (v_2(a^2-1)-1)   (a ungerade, H >= 0)
   v_p(a^^H)             = v_p(a) * a^^(H-1)              (p | a)
"""
import sys
from ar_lib import *

# --- Gegenprobe 1: tower_mod gegen direkte Ganzzahlen (unabhaengig: pow mit exakten Exponenten)
ok = 0
for a in range(2, 8):
    for n in range(0, 4):
        ex = 1
        for _ in range(n):
            ex = a ** ex if ex < 10 ** 6 else None
            if ex is None:
                break
        if ex is None:
            continue
        for m in list(range(1, 200)) + [10 ** 6, 2 ** 20, 3 ** 12, 7 ** 7]:
            assert tower_mod(a, n, m) == ex % m, (a, n, m)
            ok += 1
# 2^^5 = 2^65536 exakt, 3^^4 mod m ueber phi-Methode (unabhaengige Implementierung)
from sympy import totient
def tower_phi(a, n, m):
    # Gegenprobe mit phi statt lambda und Offset phi(m)
    if m == 1:
        return 0
    if n == 0:
        return 1 % m
    if n == 1:
        return a % m
    ph = int(totient(m))
    e = tower_phi(a, n - 1, ph)
    small = tower_exact(a, n - 1)
    if small is not None and small < ph:
        return pow(a, small, m)
    return pow(a, e + ph, m)
two5 = 2 ** 65536
for m in range(1, 3000):
    assert tower_mod(2, 5, m) == two5 % m
    for a in (3, 5, 6, 7, 10, 12):
        for n in (4, 5, 6):
            assert tower_mod(a, n, m) == tower_phi(a, n, m), (a, n, m)
    ok += 1
print("Gegenprobe tower_mod:", ok, "Faelle ok")

# --- Bewertungsformel
def alpha_mod(a, p, K):
    return tower_mod(a, 4 * K + 40, p ** K)

viol = 0; checked = 0
K = 40
for p in (2, 3, 5, 7, 11, 13):
    for a in range(2, 400):
        al = alpha_mod(a, p, K)
        for H in range(0, 30):
            x = tower_mod(a, H, p ** K)
            obs = vp((x - al) % p ** K, p)
            if a % p == 0:
                ex = tower_exact(a, H - 1)
                pred = vp(a, p) * ex if ex is not None else 10 ** 9
            elif p == 2:
                pred = vp(a - 1, 2) + H * (vp(a * a - 1, 2) - 1)
            else:
                V = vp(pow(a, p - 1) - 1, p)
                t = stab_height(a, order_mod(a, p))
                pred = max(0, V * (H - t))
                if H < t:
                    continue
            if pred >= K and obs >= K:
                continue
            checked += 1
            if obs != pred:
                viol += 1
                if viol < 10:
                    print("ABWEICHUNG", p, a, H, obs, pred)
print("Bewertungsformel: geprueft", checked, "Faelle, Abweichungen", viol)
