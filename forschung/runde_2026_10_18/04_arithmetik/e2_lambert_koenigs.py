"""E2: (a) p-adischer Turmgrenzwert = Lambert-Reihe  omega^c * sum (n+1)^{n-1} (omega^c log_p a)^n / n!
       (b) Isometrie: v_p(alpha(a') - alpha(a)) = v_p(a'-a) bei gleicher Klasse
       (c) Multiplikator lambda = log_p alpha, Koenigs-Konstante kappa = lim (a^^n - alpha)/lambda^n
"""
from ar_lib import *

K = 60
print("(a) Lambert-Reihe gegen Turm (", K, "p-adische Stellen)")
bad = 0; n_ok = 0
for p in (2, 3, 5, 7, 11, 13, 101):
    for a in range(2, 120):
        if a % p == 0:
            continue
        tow = tower_mod(a, 6 * K + 40, p ** K)
        ser = alpha_series(a, p, K)
        if tow != ser:
            bad += 1
            print("  FEHLER", p, a, vp((tow - ser) % p ** K, p))
        else:
            n_ok += 1
print("  uebereinstimmend:", n_ok, " abweichend:", bad)

# Gegenprobe mit anderer Praezision K2=25 (unabhaengige Reihenabbruchstellen)
for p, a in ((5, 2), (3, 2), (7, 10), (2, 3), (2, 7)):
    assert alpha_series(a, p, 25) == tower_mod(a, 200, p ** 25)
print("  Gegenprobe K=25 ok")

print("(b) Isometrie")
viol = 0; cnt = 0
for p in (3, 5, 7):
    period = 1
    # gleiche Klasse: a' ≡ a mod p^k * Lambda, Lambda = kgV der Carmichael-Kette von p-1 und p-1
    m = p - 1; Lam = 1
    while m > 1:
        Lam = Lam * m // gcd(Lam, m); m = carmichael(m)
    for a in range(2, 60):
        if a % p == 0:
            continue
        al = tower_mod(a, 400, p ** 30)
        for k in range(1, 12):
            for b in (1, 2):
                a2 = a + b * p ** k * Lam
                al2 = tower_mod(a2, 400, p ** 30)
                obs = vp((al2 - al) % p ** 30, p)
                cnt += 1
                if obs != vp(a2 - a, p):
                    viol += 1
print("  Faelle", cnt, "Abweichungen", viol)

print("(c) Koenigs-Konstante kappa_p(a)")
N = 120
for p, a in ((5, 2), (5, 3), (5, 7), (3, 2), (3, 5), (7, 2), (7, 3), (2, 3), (2, 5), (2, 7), (11, 2)):
    q = p ** N
    al = tower_mod(a, 10 * N, q)
    lam = padic_log(al if p != 2 else al, p, N)  # log_p alpha (Iwasawa); alpha = Einheit
    # Kontrolle: lambda = alpha * log_p a
    lam2 = al * padic_log(a, p, N) % q
    assert lam == lam2, (p, a)
    V = vp(lam, p)
    lu = lam // p ** V  # Einheitsteil
    rows = []
    prev = None
    for n in range(2, 30):
        d = (tower_mod(a, n, q) - al) % q
        v = vp(d, p)
        if v + 15 > N:
            break
        # kappa_n = d / lambda^n  -> Bewertung v - nV, Einheitsteil (d/p^v) / lu^n
        prec = N - v
        ku = (d // p ** v) * pow(pow(lu, n, p ** prec), -1, p ** prec) % p ** prec
        kv = v - n * V
        if prev is not None:
            dif = vp((ku - prev[1]) % p ** min(prec, prev[2]), p) if kv == prev[0] else 0
        else:
            dif = None
        rows.append((n, kv, dif))
        prev = (kv, ku, prec)
    print(f"  p={p} a={a} V={V} v(kappa)={rows[-1][1]}  Konvergenz v(kappa_n+1 - kappa_n):", [r[2] for r in rows[1:12]])
    print(f"      kappa (Einheitsteil, die letzten 12 {p}-adischen Ziffern mod p^12) =", prev[1] % p ** 12)
