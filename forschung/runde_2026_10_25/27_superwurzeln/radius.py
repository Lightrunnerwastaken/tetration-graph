# Konvergenzradius R_n der Taylorreihe von omega_n um 0 (Hauptzweig).
# Methode: Fortsetzung von omega_n entlang Strahlen y = s e^{i phi} (Newton auf eps_n(x) = y, Start x=0),
# Minimum von |eps_n'(x(s))| merkt sich die Naehe zu einem Verzweigungspunkt; Kandidat wird mit
# Newton auf eps_n' = 0 poliert (mpmath 30 Stellen) und dann sichtbar gemacht: Fortsetzung genau auf
# den Strahl zum kritischen Wert v = eps_n(c) muss in c enden.
import sys, cmath, math
import mpmath as mp
mp.mp.dps = 30

def ed(n, x):
    e, d = x, 1.0
    for _ in range(n):
        E = cmath.exp(e)
        e, d = x * E, E + x * E * d
    return e, d

def ed3(n, x):
    e, d, dd = x, mp.mpf(1), mp.mpf(0)
    for _ in range(n):
        E = mp.exp(e)
        e, d, dd = x * E, E + x * E * d, 2 * E * d + x * E * (d * d + dd)
    return e, d, dd

def follow(n, target, steps=3000, stop_at=1.0):
    """x(s) mit eps_n(x)=s*target, s in [0, stop_at]; liefert Endpunkt und Pfadminimum |eps'|."""
    x = 0j
    best = (1e300, None, None)
    for i in range(1, steps + 1):
        s = stop_at * i / steps
        y = s * target
        for _ in range(8):
            e, d = ed(n, x)
            if d == 0:
                return None, best
            dx = (e - y) / d
            x -= dx
            if abs(dx) < 1e-15:
                break
        e, d = ed(n, x)
        if abs(d) < best[0]:
            best = (abs(d), x, s)
    return x, best

def radius(n, nphi=720, smax=1.3):
    cands = []
    for k in range(nphi):
        phi = 2 * math.pi * (k + 0.5) / nphi
        x, best = follow(n, smax * cmath.exp(1j * phi), steps=1500)
        cands.append((best[0], best[1], best[2] * smax if best[2] else None, phi))
    cands.sort(key=lambda r: r[0])
    seen = []
    out = []
    for bd, xb, s, phi in cands[:40]:
        if xb is None:
            continue
        # polieren
        z = mp.mpc(xb)
        try:
            for _ in range(60):
                e, d, dd = ed3(n, z)
                st = d / dd
                z -= st
                if abs(st) < mp.mpf(10) ** -25:
                    break
        except Exception:
            continue
        if any(abs(z - w) < 1e-10 for w in seen):
            continue
        seen.append(z)
        v = ed3(n, z)[0]
        # Sichtbarkeit: Fortsetzung bis 0.999999 v
        xe, _ = follow(n, complex(v), steps=4000, stop_at=1 - 1e-6)
        vis = xe is not None and abs(xe - complex(z)) < 1e-2
        out.append((abs(v), z, v, vis))
    out.sort(key=lambda r: r[0])
    return out

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    for n in ns:
        out = radius(n)
        vis = [r for r in out if r[3]]
        print(f"n={n}: R_n = {mp.nstr(vis[0][0], 15) if vis else '?'}")
        for r in out[:6]:
            print(f"   c={mp.nstr(r[1],12)}  v={mp.nstr(r[2],12)}  |v|={mp.nstr(r[0],12)}  sichtbar={r[3]}")
        sys.stdout.flush()
