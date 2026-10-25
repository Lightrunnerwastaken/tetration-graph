# Kritische Punkte von eps_n (= Verzweigungspunkte von omega_n) und kritische Werte.
# eps_n'(x) = exp(eps_{n-1}) * (1 + x eps_{n-1}'(x))  =>  kritisch  <=>  1 + x eps_{n-1}'(x) = 0.
import sys
import mpmath as mp
mp.mp.dps = 30

def eps_d(n, x):
    """eps_n(x), eps_n'(x), eps_n''(x)"""
    e, d, dd = x, mp.mpf(1), mp.mpf(0)
    for _ in range(n):
        if mp.re(e) > 700:
            raise OverflowError
        E = mp.exp(e)
        # eps_{m+1} = x E, E=exp(eps_m)
        e1 = x * E
        d1 = E + x * E * d
        dd1 = 2 * E * d + x * E * (d * d + dd)
        e, d, dd = e1, d1, dd1
    return e, d, dd

def crit(n, R=2.5, grid=41):
    found = []
    for i in range(grid):
        for j in range(grid):
            z0 = mp.mpc(-R + 2 * R * i / (grid - 1), -R + 2 * R * j / (grid - 1))
            z = z0
            ok = False
            for it in range(60):
                try:
                    e, d, dd = eps_d(n, z)
                except OverflowError:
                    break
                if abs(dd) == 0:
                    break
                step = d / dd
                z = z - step
                if abs(z) > 6:
                    break
                if abs(step) < mp.mpf(10) ** -25:
                    ok = True
                    break
            if ok and abs(z) < R * 1.2:
                if all(abs(z - w) > 1e-12 for w in found):
                    found.append(z)
    return found

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(1, NMAX + 1):
        cs = crit(n)
        rows = []
        for c in cs:
            v = eps_d(n, c)[0]
            rows.append((abs(v), c, v))
        rows.sort(key=lambda r: r[0])
        print(f"n={n}: {len(cs)} kritische Punkte in |x|<~3; kleinste |kritischer Wert|:")
        for r in rows[:6]:
            print(f"   c={mp.nstr(r[1],12)}  eps_n(c)={mp.nstr(r[2],12)}  |.|={mp.nstr(r[0],12)}")
