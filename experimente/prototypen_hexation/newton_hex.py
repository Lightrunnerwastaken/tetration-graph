"""Kneser-artige Hexation zur Basis 2: Newton für H₀(y) = H(iy) im Rechteck [0, 2] × [−A, A], nur Vorwärts-Pentation mit Blattverfolgung.

 Mittellinie: H(1 + iy) = P(H₀(y)) muss gleich dem Cauchy-Integral über den Rand sein:
   links H₀, rechts H₂ = P(P(H₀)), oben q, unten q̄  (Symmetrie H(z̄) = conj H(z)).
 Normierung H(0) = 1 als Zusatzgleichung. P entlang des Weges y ↦ H₀(y) (Start y = 0, reell) mit Blattverfolgung (blatt.py);
 P′ aus Differenzen entlang des Weges. Aufruf: newton_hex.py A N dy iter
"""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
from blatt import BlattPentation, verfolge
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

A, N, DY, ITER = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
BP = BlattPentation("2")
q = complex(mp.mpc("-1.992948435932322332933652", "1.445128496639205463533289"))

xg, wg = np.polynomial.legendre.leggauss(N)
yk = A * xg                    # alle Knoten (symmetrisch)
wk = A * wg
pos = np.where(yk > 0)[0]
yp = yk[pos]
m = len(yp)
neg_zu_pos = np.array([pos[np.argmin(np.abs(yk[pos] + y))] for y in yk])   # Index des Spiegelknotens

# feines Gitter für die Wege, Knoten eingefügt
fein = np.unique(np.concatenate([np.arange(0, A + 1e-12, DY), yp, [A]]))
idx_knoten = np.array([np.argmin(np.abs(fein - y)) for y in yp])


def interp_knoten(werte_pos, ziel):
    """Polynominterpolation (baryzentrisch, Legendre-Knoten) von H₀ auf ziel; werte_pos an yp, Spiegel konjugiert."""
    alle = np.empty(N, complex)
    alle[pos] = werte_pos
    neg = np.where(yk < 0)[0]
    alle[neg] = np.conj(werte_pos[np.searchsorted(yp, -yk[neg])])
    # baryzentrische Gewichte für Legendre-Knoten
    lam = np.array([1.0 / np.prod([(yk[i] - yk[j]) / A for j in range(N) if j != i]) for i in range(N)])
    out = np.empty(len(ziel), complex)
    for t, z in enumerate(ziel):
        d = z - yk
        if np.min(np.abs(d)) < 1e-14:
            out[t] = alle[np.argmin(np.abs(d))]
        else:
            r = lam / d
            out[t] = (r @ alle) / r.sum()
    return out


def deckel(z):
    return q / (2j * np.pi) * (np.log(0 + 1j * A - z) - np.log(2 + 1j * A - z)) + np.conj(q) / (2j * np.pi) * (np.log(2 - 1j * A - z) - np.log(0 - 1j * A - z))


zM = 1 + 1j * yp                     # Auswertepunkte auf der Mittellinie
D = deckel(zM)
KR = wk[None, :] / (1 + 1j * (yk[None, :] - yp[:, None])) / (2 * np.pi)      # rechter Rand (Re = 2) relativ zu Re = 1
KL = wk[None, :] / (-1 + 1j * (yk[None, :] - yp[:, None])) / (2 * np.pi)     # linker Rand (Re = 0)


def ableitung_weg(werte, y):
    """dW/dy auf dem feinen Gitter (zentral, 5 Punkte, ungleichmässig → numpy.gradient zweiter Ordnung als Rückfall)."""
    return np.gradient(werte, y, edge_order=2)


SPRUNG = [0.0]


def auswerten(H0p):
    """H₁, H₂ und Ableitungen P′(H₀), P′(H₁) an den positiven Knoten (Weg y ↦ H₀(y) adaptiv verfeinert)."""
    def stufen(param):
        y = np.array([float(t) for t in param])
        h0 = interp_knoten(H0p, y)
        h0[0] = h0[0].real                                 # Start reell
        h1, s1 = BP.P_pfad([mp.mpc(h) for h in h0])
        h2, s2 = BP.P_pfad(h1)
        return list(zip(h0, h1, h2)), [max(a, b) for a, b in zip(s1, s2)]
    param, werte, smax = verfolge(list(fein), stufen, log=lambda *a: None)
    SPRUNG[0] = float(smax)
    y = np.array([float(t) for t in param])
    H0f = np.array([complex(w[0]) for w in werte]); H1f = np.array([complex(w[1]) for w in werte]); H2f = np.array([complex(w[2]) for w in werte])
    d1 = ableitung_weg(H1f, y) / ableitung_weg(H0f, y)
    d2 = ableitung_weg(H2f, y) / ableitung_weg(H1f, y)
    ik = np.array([np.argmin(np.abs(y - t)) for t in yp])
    return H0f, H1f, H2f, d1[ik], d2[ik], ik


def alle_knoten(werte_pos):
    out = np.empty(N, complex)
    out[pos] = werte_pos
    neg = np.where(yk < 0)[0]
    out[neg] = np.conj(werte_pos[np.searchsorted(yp, -yk[neg])])
    return out


# Startwert
st = json.load(open(HIER + f"/prototypen_hexation/start_verfolgt_{int(A)}.json"))
sy = np.array(st["y"]); sH = np.array([complex(a, b) for a, b in st["H"]])
H0 = np.interp(yp, sy, sH.real) + 1j * np.interp(yp, sy, sH.imag)
# Gewichte für H₀(0)
e0 = np.zeros(m, complex)
for j in range(m):
    v = np.zeros(m, complex); v[j] = 1.0
    e0[j] = interp_knoten(v, np.array([0.0]))[0]          # = 2 ℓ_j(0) (reell)
for it in range(ITER):
    H0f, H1f, H2f, d1, d2, ik = auswerten(H0)
    H1 = H1f[ik]; H2 = H2f[ik]
    R = H1 - (KR @ alle_knoten(H2) - KL @ alle_knoten(H0) + D)
    norm0 = (e0 @ H0).real - 1
    log(f"Iteration {it}: max |R| = {np.max(np.abs(R)):.3e}, H₀(0) − 1 = {norm0:.2e}, grösster Winkelsprung {SPRUNG[0]:.3f}, "
        f"P(H₀) oben − q = {abs(H1[-1] - q):.2e}")
    if np.max(np.abs(R)) < 1e-12 and abs(norm0) < 1e-12:
        break
    # Jacobi-Matrix (reell, 2m × 2m) plus Normierungszeile
    J = np.zeros((2 * m + 1, 2 * m))
    for j in range(m):
        for kind, e in ((0, 1.0), (1, 1j)):
            dH0 = np.zeros(m, complex); dH0[j] = e
            dH1 = d1 * dH0
            dH2 = d2 * d1 * dH0
            dR = dH1 - (KR @ alle_knoten(dH2) - KL @ alle_knoten(dH0))
            J[:m, 2 * j + kind] = dR.real
            J[m:2 * m, 2 * j + kind] = dR.imag
            J[2 * m, 2 * j + kind] = (e0 @ dH0).real
    rhs = -np.concatenate([R.real, R.imag, [norm0]])
    dx, *_ = np.linalg.lstsq(J, rhs, rcond=None)
    schritt = dx[0::2] + 1j * dx[1::2]
    daempf = min(1.0, 0.5 / max(1e-30, np.max(np.abs(schritt))))
    H0 = H0 + daempf * schritt
    log(f"   Schritt max {np.max(np.abs(schritt)):.2e}, Dämpfung {daempf:.2f}")

# H im Inneren des Rechtecks [0, 2] × [−A, A] per Cauchy-Formel
H0all, H2all = alle_knoten(H0), alle_knoten(H2)
def H_innen(z):
    rechts = (wk * H2all / (2 + 1j * yk - z)).sum() / (2 * np.pi)
    links = (wk * H0all / (0 + 1j * yk - z)).sum() / (2 * np.pi)
    top = q / (2j * np.pi) * (np.log(0 + 1j * A - z) - np.log(2 + 1j * A - z)) + np.conj(q) / (2j * np.pi) * (np.log(2 - 1j * A - z) - np.log(0 - 1j * A - z))
    return rechts - links + top
werte = {x: H_innen(x) for x in (0.5, 1.0, 1.5)}
log("H(0.5), H(1), H(1.5):", {x: f"{v.real:.12f} {v.imag:+.1e}i" for x, v in werte.items()}, " (H(1) soll 2; Levenstein 2↑↑↑↑1.5 ≈ 2.6729)")
json.dump({"A": A, "N": N, "DY": DY, "H0": [[h.real, h.imag] for h in H0], "yp": list(yp),
           "werte": {str(x): [v.real, v.imag] for x, v in werte.items()}}, open(HIER + f"/prototypen_hexation/newton_hex_{int(A)}_{N}.json", "w"), indent=1)
log("fertig")
