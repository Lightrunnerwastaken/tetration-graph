"""E4: (1) Periodenpaare nahe der reellen Achse; (2) Koebe-Radius und de-Branges-Probe fuer D(0,2);
(3) Schlichtheit von tau = z + theta(z) auf {Im z > y0} (periodische Randkurve, Translate +-1, +-2)."""
import json
import numpy as np
import lib21 as S

P = S.P_PER
erg = {}


def dtau(z, h=1e-7):
    return (S.tau(z + h) - S.tau(z - h)) / (2 * h)


# (1)
tab = []
for y in [0.2, 0.15, 0.1, 0.07, 0.05, 0.035, 0.025]:
    best = None
    for x in np.linspace(0.6, 0.95, 36):
        z = complex(x, y)
        ziel = S.tau(np.array([z]))[0] + P
        zp = z + P
        for _ in range(80):
            f = S.tau(np.array([zp]))[0] - ziel
            zp = zp - f / dtau(np.array([zp]))[0]
            if abs(f) < 1e-12:
                break
        else:
            continue
        if zp.imag < 0.2:
            continue
        fe = abs(S.sexp(np.array([zp]))[0] - S.sexp(np.array([z]))[0])
        if best is None or (zp - z).real < best[0]:
            best = ((zp - z).real, x, zp, fe)
    tab.append((y, best[0], best[1], str(best[2]), best[3]))
    print(f"y={y}: min Re(z'-z)={best[0]:.6f} bei x={best[1]:.3f}, z'={best[2]:.4f}, |sexp z'-sexp z|={best[3]:.1e}", flush=True)
erg["paare_klein_y"] = tab

# (2) Koebe / Landau fuer D(0,2): Abstand von sexp(0)=1 zum Bildrand
t = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
for R in (1.99, 1.999):
    w = S.sexp(R * np.exp(1j * t))
    d = np.abs(w - 1).min()
    k = int(np.argmin(np.abs(w - 1)))
    erg[f"koebe_R{R}"] = (d, d / (R * S.TAYLOR[1]), str(R * np.exp(1j * t[k])))
    print(f"R={R}: dist(1, Rand) = {d:.6f}, Verhaeltnis zu R*sexp'(0) = {d/(R*S.TAYLOR[1]):.6f} (Koebe: >= 0.25), bei z={R*np.exp(1j*t[k]):.4f}")
c = S.TAYLOR
q = [abs(c[n]) * 2 ** (n - 1) / c[1] / n for n in range(2, 128)]
erg["debranges_max"] = max(q)
print("max_n |c_n| 2^(n-1)/(n c_1) =", max(q), "bei n =", 2 + int(np.argmax(q)), "; n|c_n|2^n fuer n=20,60,120:",
      [n * abs(c[n]) * 2 ** n for n in (20, 60, 120)])

# (3) tau-Schlichtheit
for y0 in [0.3, 0.1, 0.05, 0.03]:
    x = np.linspace(0, 1, 3000)
    g = S.tau(x + 1j * y0)
    # Kurve ueber [-2, 3] zusammensetzen (tau(z+1) = tau(z) + 1)
    kurve = np.concatenate([g + k for k in range(-2, 3)])
    dp = S.doppelpunkte(kurve, False, 1)
    erg[f"tau_y0_{y0}"] = "einfach" if not dp else str(dp)
    print(f"tau auf Im z = {y0}: {'doppelpunktfrei' if not dp else dp}; Spannweite Im tau: {g.imag.min():.3f}..{g.imag.max():.3f}", flush=True)
json.dump(erg, open("daten/e4_rest.json", "w", encoding="utf-8"), indent=1, default=str)
