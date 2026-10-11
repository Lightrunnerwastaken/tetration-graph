"""Kneser-artige Hexation, zweiter Anlauf: Kurve γ(w) = α_q(w) (F_q(γ) = w, w reell) auf Tschebyschow-Knoten, dann Fornberg-Verfahren
für die Randzuordnung w(x): γ(w(x)) − x darf nur Frequenzen n ≥ 0 haben. Aufruf: hex_fornberg.py B L0 q0re q0im M K"""
import sys, time, json
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
import numpy as np
import mpmath as mp
from hexq import HexQ
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

B, L0, q0 = sys.argv[1], mp.mpf(sys.argv[2]), mp.mpc(sys.argv[3], sys.argv[4])
M, K = int(sys.argv[5]), int(sys.argv[6])
H = HexQ(B, L0, q0)
log("q", mp.nstr(H.q, 20), "mu", mp.nstr(H.mu, 12))

# --- Startwerte aus dem groben Lauf (Kurve als Graph g(t), K=128)
alt = json.load(open(r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad/hex_lauf_" + B + ".json"))
galt = [mp.mpf(x) for x in alt["128"]["g"]]
talt = [mp.mpf(j) / 128 for j in range(128)]
Falt = [mp.re(v) for v in H.F([t + 1j * g for t, g in zip(talt, galt)])]
log("alte Kurve: F von", mp.nstr(Falt[0], 10), "bis", mp.nstr(Falt[-1], 10), "monoton:", all(Falt[i] < Falt[i + 1] for i in range(127)))
wa = Falt[0]
wb = H.P([wa])[0]
log("Grundintervall w ∈ [", mp.nstr(wa, 12), ",", mp.nstr(wb, 12), "]")
# Tschebyschow-Knoten (zweiter Art, mit Endpunkten) im w-Intervall
wk = [(wa + wb) / 2 - (wb - wa) / 2 * mp.cos(mp.pi * k / (M - 1)) for k in range(M)]


def start(w):
    # lineare Interpolation in (F_alt, t, g), periodisch fortgesetzt über P
    if w > Falt[-1]:
        return talt[-1] + 1j * galt[-1]
    i = max(i for i in range(128) if Falt[i] <= w) if w >= Falt[0] else 0
    j = min(i + 1, 127)
    s = (w - Falt[i]) / (Falt[j] - Falt[i]) if j != i else 0
    return talt[i] + s * (talt[j] - talt[i]) + 1j * (galt[i] + s * (galt[j] - galt[i]))


zk = [start(w) for w in wk]
zk[-1] = zk[0] + 1                                  # Endpunkt: γ(P(wa)) = γ(wa) + 1
h = mp.mpf("1e-12")
for it in range(12):
    v = H.F(zk + [z + h for z in zk] + [z - h for z in zk])
    f, fp_, fm = v[:M], v[M:2 * M], v[2 * M:]
    neu = [z - (fz - w) / ((a - b) / (2 * h)) for z, fz, w, a, b in zip(zk, f, wk, fp_, fm)]
    d = max(abs(a - b) for a, b in zip(neu, zk))
    zk = neu
    log(f"Newton {it}: max Schritt {mp.nstr(d, 3)}, max |F − w| {mp.nstr(max(abs(a - b) for a, b in zip(f, wk)), 3)}")
    if d < mp.mpf(10) ** -25:
        break
log("Periodizität γ(P(wa)) − γ(wa) − 1 =", mp.nstr(zk[-1] - zk[0] - 1, 3))
log("Kurve γ(w):", [mp.nstr(z, 8) for z in zk[::max(1, M // 8)]])

# --- Tschebyschow-Interpolation von γ(w) in numpy (komplex)
wk_n = np.array([float(w) for w in wk])
gk_n = np.array([complex(z) for z in zk])
a_, b_ = float(wa), float(wb)
xk = -np.cos(np.pi * np.arange(M) / (M - 1))      # in [-1, 1], gleiche Reihenfolge wie wk
cheb = np.polynomial.chebyshev.chebfit(xk, gk_n, M - 1)
log("Tschebyschow-Koeffizienten |c_n| (n = 0, M/4, M/2, 3M/4, M−1):", [f"{abs(cheb[n]):.2e}" for n in (0, M // 4, M // 2, 3 * M // 4, M - 1)])
dcheb = np.polynomial.chebyshev.chebder(cheb)


def gam(w):
    x = (2 * w - a_ - b_) / (b_ - a_)
    return np.polynomial.chebyshev.chebval(x, cheb)


def dgam(w):
    x = (2 * w - a_ - b_) / (b_ - a_)
    return np.polynomial.chebyshev.chebval(x, dcheb) * 2 / (b_ - a_)


# --- Fornberg: Unbekannte w_j = w(x_j), x_j = j/K, w_0 = wa fest
xs = np.arange(K) / K
# Startzuordnung: Re γ(w) ≈ x  (aus dem groben Lauf: t ≈ x)
wgrid = np.linspace(a_, b_, 4001)
regam = gam(wgrid).real
w = np.interp(xs, regam - regam[0], wgrid) if np.all(np.diff(regam) > 0) else a_ + (b_ - a_) * xs
w[0] = a_
n_neg = np.arange(1, K // 2)
E = np.exp(2j * np.pi * np.outer(n_neg, xs))       # e^{+2πi n x}: ĥ_{−n} = (1/K) Σ h_j e^{+2πi n x_j}
for it in range(60):
    hj = gam(w) - xs
    R = E @ hj / K                                   # Koeffizienten der negativen Frequenzen
    J = E * dgam(w)[None, :] / K
    Jr = np.vstack([J.real, J.imag])[:, 1:]
    Rr = np.concatenate([R.real, R.imag])
    dw, *_ = np.linalg.lstsq(Jr, -Rr, rcond=None)
    w[1:] += dw
    if it % 5 == 0 or np.max(np.abs(dw)) < 1e-14:
        log(f"Fornberg {it}: |R| = {np.max(np.abs(R)):.2e}, |Δw| = {np.max(np.abs(dw)):.2e}, monoton: {bool(np.all(np.diff(w) > 0))}")
    if np.max(np.abs(dw)) < 1e-14:
        break
hj = gam(w) - xs
ahat = np.fft.fft(hj) / K                            # ĥ_n, n = 0..K−1 (n ≥ K/2 sind negative Frequenzen)
log("negative Frequenzen max:", f"{np.max(np.abs(ahat[K // 2 + 1:])):.2e}", " |â_n| n=1,K/8,K/4,K/2−1:", [f"{abs(ahat[n]):.2e}" for n in (1, K // 8, K // 4, K // 2 - 1)])
apos = [mp.mpc(complex(ahat[n])) for n in range(K // 2)]


def theta(z):
    return mp.fsum(apos[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(apos)))


def Hx(xx):
    return H.F([x + theta(x) for x in xx])


probe = Hx([mp.mpf(0.5) / K, mp.mpf("0.4") + mp.mpf(1) / (3 * K), mp.mpf("0.77") + mp.mpf(1) / (5 * K)])
log("Im H zwischen den Stützstellen:", [mp.nstr(mp.im(v), 3) for v in probe])
# Normierung H(x*) = 1
xa, xb = mp.mpf("2.1"), mp.mpf("2.2")
fa, fb = mp.re(Hx([xa])[0]) - 1, mp.re(Hx([xb])[0]) - 1
for it in range(40):
    xc = xb - fb * (xb - xa) / (fb - fa)
    xa, fa, xb = xb, fb, xc
    fb = mp.re(Hx([xb])[0]) - 1
    if abs(xb - xa) < mp.mpf(10) ** -18:
        break
xs0 = xb
stellen = [mp.mpf(-1), mp.mpf("-0.5"), mp.mpf("0.5"), mp.mpf(1), mp.mpf("1.5")]
res = Hx([xs0 + s for s in stellen])
log(f"M={M} K={K}: x* = {mp.nstr(xs0, 15)};", {mp.nstr(s, 2): mp.nstr(r, 12) for s, r in zip(stellen, res)})
json.dump({"M": M, "K": K, "x_stern": mp.nstr(xs0, 20), "werte": {mp.nstr(s, 2): [mp.nstr(mp.re(r), 18), mp.nstr(mp.im(r), 4)] for s, r in zip(stellen, res)},
           "gamma_w": [mp.nstr(w_, 20) for w_ in wk], "gamma": [[mp.nstr(mp.re(z), 20), mp.nstr(mp.im(z), 20)] for z in zk],
           "w_x": [float(x) for x in w], "ahat": [[float(a.real), float(a.imag)] for a in ahat[:K // 2]]},
          open(sys.argv[0].replace(".py", f"_{B}_{M}_{K}.json"), "w"), indent=1)
log("fertig")
