"""E3: Periodenpaare. Zu z in H loese tau(z') = tau(z) + P (Newton, Start z+P), P = 2 pi i / L.
Dann ist sexp(z') = sexp(z) (F ist P-periodisch). Gemessen: Breite Re(z'-z) und Gegenprobe |sexp z' - sexp z|.
Ausserdem: Darboux-Test fuer Rechtecke {b-w < Re z < b, |Im z| < Y} ohne Schnitt.
"""
import json
import numpy as np
import lib21 as S

P = S.P_PER
erg = {"P": [P.real, P.imag]}


def dtau(z, h=1e-6):
    return (S.tau(z + h) - S.tau(z - h)) / (2 * h)


paare = []
best = (9e9, None)
for y in [0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]:
    for x in np.linspace(0, 1, 41)[:-1]:
        z = complex(x, y)
        ziel = S.tau(np.array([z]))[0] + P
        zp = z + P
        ok = False
        for _ in range(60):
            f = S.tau(np.array([zp]))[0] - ziel
            zp = zp - f / dtau(np.array([zp]))[0]
            if abs(f) < 1e-13:
                ok = True
                break
        if not ok or zp.imag < 0.2:
            continue
        d = zp - z
        fehler = abs(S.sexp(np.array([zp]))[0] - S.sexp(np.array([z]))[0])
        paare.append((x, y, d.real, d.imag, fehler))
        if d.real < best[0]:
            best = (d.real, (x, y, d.imag, fehler))
    zeile = [p for p in paare if p[1] == y]
    if zeile:
        m = min(zeile, key=lambda p: p[2])
        print(f"y={y}: min Re(z'-z) = {m[2]:.6f} bei x={m[0]:.3f} (Im diff {m[3]:.4f}), max Re = {max(p[2] for p in zeile):.6f}, max|sexp z'-sexp z| = {max(p[4] for p in zeile):.1e}", flush=True)
erg["paare"] = paare
erg["min_breite"] = best
print("kleinste Breite:", best)


# Darboux-Test fuer Rechtecke ohne Schnitt
def rand(b, w, Y, n=1500):
    a = b - w
    yy = np.concatenate([np.linspace(-Y, -1, n // 3, endpoint=False), np.linspace(-1, 1, n, endpoint=False),
                         np.linspace(1, Y, n // 3, endpoint=False)])
    rechts = b + 1j * yy
    oben = np.linspace(b, a, n // 3, endpoint=False) + 1j * Y
    links = a + 1j * yy[::-1]
    unten = np.linspace(a, b, n // 3, endpoint=False) - 1j * Y
    return np.concatenate([rechts, oben, links, unten])


tests = [(2.4, 1.0), (2.4, 2.0), (2.4, 3.0), (2.4, 4.3), (2.5, 1.0), (2.53, 1.0), (2.6, 1.0), (2.0, 3.95), (1.0, 2.9)]
dt = {}
for b, w in tests:
    if b - w <= -2:
        continue
    zz = rand(b, w, 8.0)
    ww = S.sexp(zz)
    dp = S.doppelpunkte(ww, True, 1)
    dt[f"b={b},w={w}"] = "einfach" if not dp else f"Doppelpunkt bei z={zz[dp[0][0]]:.4f} / {zz[dp[0][1]]:.4f}"
    print(b, w, dt[f"b={b},w={w}"], flush=True)
erg["darboux"] = dt
json.dump(erg, open("daten/e3_periodenpaare.json", "w", encoding="utf-8"), indent=1, default=str)
