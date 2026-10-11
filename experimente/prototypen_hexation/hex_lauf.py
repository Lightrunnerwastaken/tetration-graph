"""Kneser-artige Hexation: Kurve C mit K Stützstellen (Startwerte aus gröberer Rechnung), Theodorsen, Normierung, Werte. Aufruf: hex_lauf.py B L0 q0re q0im K1 K2 …"""
import sys, time, json
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
import mpmath as mp
from hexq import HexQ
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

B, L0, q0 = sys.argv[1], mp.mpf(sys.argv[2]), mp.mpc(sys.argv[3], sys.argv[4])
Ks = [int(k) for k in sys.argv[5:]]
H = HexQ(B, L0, q0)
log("q", mp.nstr(H.q, 25), "mu", mp.nstr(H.mu, 15))


def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]


def fourier_auswerten(fh, x):
    M = len(fh)
    s = fh[0]
    for n in range(1, (M + 1) // 2):
        s += 2 * mp.re(fh[n] * mp.expjpi(2 * n * x))
    if M % 2 == 0:
        s += mp.re(fh[M // 2] * mp.cos(mp.pi * M * x))
    return mp.re(s)


gh_alt = None
ergebnisse = {}
for K in Ks:
    ts = [mp.mpf(j) / K for j in range(K)]
    if gh_alt is None:
        g0, g1 = [mp.mpf("-1.45")] * K, [mp.mpf("-1.60")] * K
    else:
        g1 = [fourier_auswerten(gh_alt, t) for t in ts]
        g0 = [g - mp.mpf("1e-4") for g in g1]
    f0 = [mp.im(v) for v in H.F([t + 1j * g for t, g in zip(ts, g0)])]
    vals = H.F([t + 1j * g for t, g in zip(ts, g1)])
    f1 = [mp.im(v) for v in vals]
    for it in range(40):
        g2 = [b if fb == fa else b - fb * (b - a) / (fb - fa) for a, b, fa, fb in zip(g0, g1, f0, f1)]
        g0, f0, g1 = g1, f1, g2
        vals = H.F([t + 1j * g for t, g in zip(ts, g1)])
        f1 = [mp.im(v) for v in vals]
        d = max(abs(a - b) for a, b in zip(g0, g1))
        log(f"K={K} Sekante {it}: max Δg {mp.nstr(d, 3)}, max |Im F| {mp.nstr(max(abs(x) for x in f1), 3)}")
        if d < mp.mpf(10) ** -22:
            break
    g = g1
    gh = dft(g)
    gh_alt = gh
    log(f"K={K} |ĝ_n|:", [(n, mp.nstr(abs(gh[n]), 3)) for n in (1, 4, 8, 16, 32, 48, 64) if n < K // 2])
    j_min = min(range(K), key=lambda j: g[j])
    log(f"K={K} tiefster Punkt von C: t = {mp.nstr(ts[j_min], 4)}, g = {mp.nstr(g[j_min], 10)}")
    # Theodorsen
    xs = ts
    u = [mp.mpf(0)] * K
    for it in range(400):
        v = [fourier_auswerten(gh, x + uu) for x, uu in zip(xs, u)]
        vh = dft(v)
        a = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
        u_neu = [mp.re(mp.fsum(a[n] * mp.expjpi(2 * n * x) for n in range(len(a)))) for x in xs]
        dd = max(abs(p - q) for p, q in zip(u, u_neu))
        u = u_neu
        if dd < mp.mpf(10) ** -22:
            break
    log(f"K={K} Theodorsen {it} Schritte, |a_n|:", [(n, mp.nstr(abs(a[n]), 3)) for n in (1, 4, 8, 16, 32, 48, 63) if n < len(a)])

    def theta(z):
        return mp.fsum(a[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(a)))

    def Hx(xx):
        return H.F([x + theta(x) for x in xx])

    probe = Hx([mp.mpf("0.1") + mp.mpf(1) / (3 * K), mp.mpf("0.4") + mp.mpf(1) / (3 * K), mp.mpf("0.83") + mp.mpf(1) / (3 * K)])
    log(f"K={K} Im H zwischen Stützstellen:", [mp.nstr(mp.im(v), 3) for v in probe])
    xa, xb = mp.mpf("2.1"), mp.mpf("2.2")
    fa, fb = mp.re(Hx([xa])[0]) - 1, mp.re(Hx([xb])[0]) - 1
    for it in range(40):
        xc = xb - fb * (xb - xa) / (fb - fa)
        xa, fa, xb = xb, fb, xc
        fb = mp.re(Hx([xb])[0]) - 1
        if abs(xb - xa) < mp.mpf(10) ** -20:
            break
    xs0 = xb
    stellen = [mp.mpf(-1), mp.mpf("-0.5"), mp.mpf("0.5"), mp.mpf(1), mp.mpf("1.5")]
    res = Hx([xs0 + s for s in stellen])
    log(f"K={K} x* = {mp.nstr(xs0, 15)};", {mp.nstr(s, 2): mp.nstr(r, 12) for s, r in zip(stellen, res)})
    ergebnisse[K] = {"x_stern": mp.nstr(xs0, 25), "werte": {mp.nstr(s, 2): [mp.nstr(mp.re(r), 20), mp.nstr(mp.im(r), 5)] for s, r in zip(stellen, res)},
                     "g": [mp.nstr(x, 25) for x in g]}
    json.dump(ergebnisse, open(sys.argv[0].replace(".py", f"_{B}.json"), "w"), indent=1)
log("fertig")
