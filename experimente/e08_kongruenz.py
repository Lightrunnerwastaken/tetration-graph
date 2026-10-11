"""exp.kongruenz — stabile Endziffern von a↑↑n und Kongruenzgeschwindigkeit V(a).

Prüft br.satz.turm_mod_m (Konstanz mod m) und die zitierte Formel V(a) = v_2(a²−1) − 1 für a ≡ 5 (mod 10)
(OEIS A317905) als Gegenbeispielsuche für a ≤ 400.
"""
from math import gcd

from gemeinsam import Lauf, plt

lauf = Lauf("exp.kongruenz", __file__)
K = 40               # Stellen
M = 10 ** K


def carmichael(m):
    def lam_pp(p, e):
        if p == 2:
            return 1 if e == 1 else (2 if e == 2 else 2 ** (e - 2))
        return (p - 1) * p ** (e - 1)
    res, n, p = 1, m, 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            l = lam_pp(p, e)
            res = res * l // gcd(res, l)
        p += 1
    if n > 1:
        l = n - 1
        res = res * l // gcd(res, l)
    return res


def turm_echt(a, n, grenze=10 ** 6):
    """a↑↑n exakt, falls < grenze, sonst None."""
    x = 1
    for _ in range(n):
        if x > 64:
            return None
        x = a ** x
        if x >= grenze:
            return None
    return x


def turm_mod(a, n, m):
    if m == 1:
        return 0
    if n == 0:
        return 1 % m
    lam = carmichael(m)
    e_echt = turm_echt(a, n - 1, grenze=10 ** 4)
    if e_echt is not None:
        return pow(a, e_echt, m)
    e = turm_mod(a, n - 1, lam)
    return pow(a, e + lam, m)   # Exponent gross: a^E ≡ a^{(E mod λ) + λ} (mod m) für E ≥ log2 m


def stabile(a, n):
    x, y = turm_mod(a, n, M), turm_mod(a, n + 1, M)
    s = 0
    while s < K and (x // 10 ** s) % 10 == (y // 10 ** s) % 10:
        s += 1
    return s


def v(p, x):
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


# Konsistenz der Modulrechnung gegen exakte Türme
ok = all(turm_mod(a, n, 10 ** 6) == turm_echt(a, n, 10 ** 30) % 10 ** 6 for a in range(2, 8) for n in range(1, 4)
         if turm_echt(a, n, 10 ** 30) is not None)
lauf.pruefe("Modulrechnung = exakte Türme (kleine Fälle)", ok)

geschw = {}
for a in range(2, 401):
    if a % 10 == 0:
        continue
    s = [stabile(a, n) for n in range(2, 9)]
    d = [s[i + 1] - s[i] for i in range(len(s) - 1) if s[i + 1] < K]
    if len(d) >= 2 and d[-1] == d[-2]:
        geschw[a] = d[-1]
lauf.pruefe("Kongruenzgeschwindigkeit stabilisiert für ≥ 95 % der a ≤ 400", len(geschw) >= 0.95 * 360, len(geschw))
fuenf = [a for a in geschw if a % 10 == 5]
abw = [a for a in fuenf if geschw[a] != v(2, a * a - 1) - 1]
lauf.pruefe(f"Formel V(a) = v₂(a²−1) − 1 für alle {len(fuenf)} Basen a ≡ 5 (mod 10)", not abw, abw[:10])
lauf.daten = {"V": geschw, "abweichungen_a_eq_5": abw}

fig, ax = plt.subplots(figsize=(7, 3.6))
for rest, farbe in [(1, 0), (3, 1), (5, 2), (7, 3), (9, 4), (2, 5), (4, 6), (6, 7)]:
    xs = [a for a in sorted(geschw) if a % 10 == rest]
    ax.plot(xs, [geschw[a] for a in xs], "o", ms=3, label=f"a ≡ {rest} (mod 10)")
ax.set_xlabel("Basis a")
ax.set_ylabel("V(a)")
ax.set_title("Kongruenzgeschwindigkeit: neue stabile Endziffern pro Turmstufe")
ax.legend(fontsize=7, ncol=4)
lauf.bild(fig)
lauf.zusammenfassung = (f"V(a) für {len(geschw)} Basen bestimmt (40 Stellen). Formel v₂(a²−1)−1 für a ≡ 5 (mod 10): "
                        f"{len(fuenf) - len(abw)}/{len(fuenf)} Treffer.")
lauf.schreibe()
