"""E4: exponentielle Fakultaet EF(n) = n^EF(n-1), EF(1) = 1, modulo m.
 (a) EF(n) mod m ist schliesslich periodisch; kleinste Periode P(m) gegen Lambda(m) = kgV(m, lambda(m), lambda^2(m), ...)
 (b) Vorperiode N(m) gegen Kettenlaenge l(m)
 (c) profinite Fortsetzung E^(x): Werte bei x = 0, 1, 2, 3, -1 und Vergleich mit EF(n)
"""
from functools import lru_cache
from ar_lib import carmichael, nu, vp
from math import gcd

EF_small = {0: 1, 1: 1, 2: 2, 3: 9, 4: 4 ** 9}  # EF(0) := 1 (leeres Produkt-Konvention, nur Hilfswert)


@lru_cache(maxsize=None)
def EF_mod(n, m):
    if m == 1:
        return 0
    if n <= 4:
        return EF_small[n] % m
    # Exponent EF(n-1) >= EF(4) = 262144 >= nu(m) fuer alle hier benutzten m
    lam = carmichael(m)
    e = EF_mod(n - 1, lam)
    return pow(n, e % lam + lam * (nu(m) + 1), m)


def Lam(m):
    L = 1
    while m > 1:
        L = L * m // gcd(L, m)
        m = carmichael(m)
    return L


def chainlen(m):
    l = 0
    while m > 1:
        m = carmichael(m); l += 1
    return l


def main():
    # Gegenprobe EF_mod gegen exakte Werte EF(5) = 5^262144 und EF(6) mod m ueber phi-unabhaengig
    EF5 = 5 ** 262144
    for m in range(1, 3000):
        assert EF_mod(5, m) == EF5 % m
    print("Gegenprobe EF(5) ok")

    res = []
    mismatch = []
    skipped = []
    for m in range(2, 401):
        L = Lam(m)
        if L > 60000:
            skipped.append(m); continue
        # Folge fuer n = 5 .. 5 + 3L + 40 ; Periodensuche
        top = 6 + 3 * L + 40
        seq = [EF_mod(n, m) for n in range(5, top)]
        # kleinste Periode P (Teiler von L) mit seq[i] == seq[i+P] fuer i >= i0 (i0 klein)
        best = None
        divs = [d for d in range(1, L + 1) if L % d == 0]
        for P in divs:
            # Vorperiode: kleinstes i0 mit Periodizitaet ab i0 (in der ganzen Folge)
            i0 = None
            for start in range(0, 40):
                if all(seq[i] == seq[i + P] for i in range(start, len(seq) - P)):
                    i0 = start; break
            if i0 is not None:
                best = (P, i0 + 5); break
        P, N0 = best
        res.append((m, P, L, N0, chainlen(m)))
        if P != L:
            mismatch.append((m, P, L))
    print("uebersprungen (Lambda>60000):", len(skipped))
    print("m in [2,400]: P(m) = Lambda(m) in", sum(1 for r in res if r[1] == r[2]), "von", len(res), "Faellen")
    print("  Ausnahmen (m, P, Lambda), erste 25:", mismatch[:25])
    # Vorperiode: EF(n) mod m periodisch ab N0(m); Vergleich mit Kettenlaenge (Vorperiode ab n=5 gezaehlt)
    mx = max(r[3] for r in res)
    print("  groesste Vorperiode N0 (ab n>=5 gemessen):", mx, "; Faelle mit N0 > 5:", [(r[0], r[3], r[4]) for r in res if r[3] > 5][:20])

    # (c) profinite Werte: E^(x) mod m = EF(n) mod m, n ≡ x mod Lambda(m), n gross
    def Ehat(x, m):
        L = Lam(m)
        n = x % L + L * 50 + 0
        while n < 60:
            n += L
        return EF_mod(n, m)

    for m in (2 ** 20, 3 ** 12, 5 ** 8, 7 ** 6, 10 ** 6):
        vals = {x: Ehat(x, m) for x in (-2, -1, 0, 1, 2, 3, 4, 5, 6)}
        ef = {x: EF_mod(x, m) for x in (1, 2, 3, 4, 5, 6)}
        print(f"  m={m}: E^(x) =", vals)
        print(f"          EF(x) =", ef)


if __name__ == "__main__":
    main()
