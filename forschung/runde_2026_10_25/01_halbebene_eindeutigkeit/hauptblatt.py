# Hauptblatt-Singularitaeten von T = sexp_K o sigma, sigma = id + p, p reelles trigonometrisches Polynom.
# Sucht z* in H = {Re z > -2}, Im z* > 0 mit sigma(z*) = -2, das von der reellen Achse aus auf einem Weg erreicht wird,
# auf dem sigma die Schnittlinie (-inf, -2] nicht trifft (Hauptblatt von sexp_K). Dann hat T bei z* eine Log-Singularitaet.
# Kein Kneser-Wert noetig: nur sigma wird ausgewertet. Aufruf: python hauptblatt.py [anzahl_zufall] [seed]
import sys, random
import mpmath as mp
mp.mp.dps = 30

def mk(coef):
    # coef: Liste (n, a_n, b_n): p(z) = sum a_n (cos(2 pi n z) - 1) + b_n sin(2 pi n z)  => p(0) = 0, reell
    def p(z):
        return sum(a*(mp.cos(2*mp.pi*n*z)-1) + b*mp.sin(2*mp.pi*n*z) for n, a, b in coef)
    def dp(z):
        return sum(2*mp.pi*n*(-a*mp.sin(2*mp.pi*n*z) + b*mp.cos(2*mp.pi*n*z)) for n, a, b in coef)
    return p, dp

def monoton(coef):
    _, dp = mk(coef)
    return min(1 + dp(mp.mpf(k)/400) for k in range(400)) > 0.05

def weg_ok(sig, z, schritte=4000):
    # senkrechter Weg von Re z* (reell) nach z*; prueft Vorzeichenwechsel von Im sigma vor dem Ende
    x, y = mp.re(z), mp.im(z)
    alt = None
    for k in range(1, schritte):
        w = sig(mp.mpc(x, y*k/schritte))
        if alt is not None and mp.sign(mp.im(w)) != mp.sign(mp.im(alt)):
            # Kreuzung der reellen Achse: nur erlaubt bei Re sigma > -2
            if mp.re(w) <= -2:
                return False
        alt = w
    return True

def suche(coef, starts=60, rnd=None):
    p, dp = mk(coef)
    sig = lambda z: z + p(z)
    dsig = lambda z: 1 + dp(z)
    rnd = rnd or random.Random(1)
    funde = []
    for _ in range(starts):
        z = mp.mpc(rnd.uniform(-1.99, 0.0), rnd.uniform(0.05, 2.5))
        try:
            for _ in range(60):
                dz = (sig(z) + 2)/dsig(z)
                z -= dz
                if abs(dz) < mp.mpf(10)**-25:
                    break
        except Exception:
            continue
        if abs(sig(z)+2) < 1e-20 and mp.re(z) > -2 and mp.im(z) > 1e-6 and abs(dsig(z)) > 1e-8:
            if all(abs(z-f) > 1e-10 for f in funde):
                funde.append(z)
    funde.sort(key=lambda z: mp.im(z))
    for z in funde:
        if weg_ok(sig, z):
            return z, abs(dsig(z)), len(funde)
    return None, None, len(funde)

if __name__ == '__main__':
    beispiele = {
        'eps sin, eps=0.1': [(1, 0, 0.1)],
        'eps sin, eps=0.01': [(1, 0, 0.01)],
        'eps (cos-1), eps=0.1': [(1, 0.1, 0)],
        'sin 4pi z, 0.03': [(2, 0, 0.03)],
        'gemischt': [(1, 0.05, 0.08), (3, -0.01, 0.005)],
    }
    for name, c in beispiele.items():
        z, d, n = suche(c)
        print(f'{name:22s} monoton={monoton(c)}  Wurzeln={n:3d}  z*={mp.nstr(z, 12) if z is not None else None}  |sigma\'(z*)|={mp.nstr(d, 6) if d else None}')
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    rnd = random.Random(seed)
    ok = keine = 0
    for i in range(N):
        grad = rnd.randint(1, 4)
        c = [(n, rnd.uniform(-1, 1)*0.15/n**1.5, rnd.uniform(-1, 1)*0.15/n**1.5) for n in range(1, grad+1)]
        if not monoton(c):
            continue
        z, d, n = suche(c, starts=40, rnd=random.Random(seed*1000+i))
        if z is None:
            keine += 1
            print('KEIN Hauptblatt-Treffer:', c, 'Wurzeln', n)
        else:
            ok += 1
    if N:
        print(f'Zufall (seed {seed}): {ok} mit Hauptblatt-Singularitaet, {keine} ohne Treffer')
