"""Gebiet 02: Schwankung von arg a_n der θ-Welle (Basis e oder 2).

Teil A  Φ(s) = arg(1+θ'(s+iη)) aus den gemessenen a_n (Etappe 2), zwei Einstellungen; Nullstellen von 1+θ' im Kreis; Mittelwerte.
Teil B  Abgleich auf der reellen Achse: s + θ(s) = h(x+i0) hat eine reelle Lösung s = slog x, und (1+θ'(s))/h'(x+i0) = sexp'(s) ist reell positiv.
Teil C  a_n gegen h'(n+i0) und gegen den führenden Term −h'(c−ln n)/n, n ≤ 4096.
Aufruf: python p1_phase.py [e|2]
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

from kern import Fix

B = sys.argv[1] if len(sys.argv) > 1 else "e"
E2 = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\identitaeten_2026_10_14\etappe2\daten")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 30
fx = Fix(mp.e if B == "e" else mp.mpf(B), K=40)
lb, lam = fx.lb, fx.lam
per = 2j * mp.pi / mp.log(lam)
kor = json.loads((E2 / f"korrektur_{B}.json").read_text())
A = mp.mpf(kor["A"])
c = mp.log(1j * A / (2 * mp.pi * lb))


def lade(fn):
    d = json.loads((E2 / fn).read_text())
    return np.array([complex(float(r), float(i)) for r, i in d["a"]])


a1 = lade(f"a_{B}_80_0.003_16384.json")
a2 = lade(f"a_{B}_60_0.002_16384.json")
N = len(a1) - 1
nn = np.arange(N + 1)


# ---------- h und h' auf der oberen Halbebene (Hauptzweige), Kettenregel ----------
def h_hp(w):
    """(h(w), h'(w)) mit h = reguläre Abelfunktion an L;  h'(w) = Π 1/(w_k ln b) · σ'(u)/(σ(u) log λ)."""
    w = mp.mpc(w)
    n, d = 0, mp.mpc(1)
    while abs(w - fx.L) > mp.mpf("1e-3") or n < 5:
        d /= w * lb
        w = mp.log(w) / lb
        n += 1
    u = w - fx.L
    tol = mp.mpf(10) ** (-mp.mp.dps)
    while abs(u) ** (fx.K + 1) > tol * abs(u):
        d /= (fx.L + u) * lb
        u = mp.log(fx.L + u) / lb - fx.L
        n += 1
    S = fx.S
    su = mp.polyval(S[::-1], u)
    dsu = mp.polyval([k * S[k] for k in range(len(S) - 1, 0, -1)], u)
    return n + mp.log(su) / mp.log(lam), d * dsu / su / mp.log(lam)


def hp_real(x):
    return h_hp(mp.mpc(x, mp.mpf("1e-22")))


# ---------- θ, θ' aus der Reihe ----------
def theta_d(a, z):
    q = np.exp(2j * np.pi * z * nn)
    return np.sum(a * q), 2j * np.pi * np.sum(nn * a * q)


def theta_d_extra(a, s, etas=(0.002, 0.003, 0.004)):
    """Werte auf der Achse durch quadratische Extrapolation η → 0 (θ ist dort analytisch, s ∉ ℤ)."""
    vs = [theta_d(a, s + 1j * e) for e in etas]
    e = np.array(etas)
    out = []
    for k in range(2):
        y = np.array([v[k] for v in vs])
        V = np.vander(e, 3)
        coef = np.linalg.solve(V, y)
        out.append(coef[-1])
    return out


res = {"basis": B}
# ================= Teil A =================
print("== Teil A: Φ(s) = arg(1+θ'(s+iη)) ==")
M = 1 << 15
for name, a, eta in [("dps80,y=.003", a1, 0.003), ("dps60,y=.002", a2, 0.002)]:
    r = np.exp(-2 * np.pi * eta)
    coef = np.zeros(M, complex)
    coef[1:N + 1] = 2j * np.pi * nn[1:] * a[1:] * r ** nn[1:]
    coef[0] = 1.0
    G = np.fft.ifft(coef) * M          # G(r e^{2πi j/M}) = 1+θ'(j/M + iη)
    ph = np.unwrap(np.angle(G))
    wind = (ph[-1] - ph[0] + np.angle(G[0] / G[-1])) / (2 * np.pi)
    print(f"{name}: Umlaufzahl von 1+θ' auf |q|=e^(-2πη): {wind:+.6f}  (0 ⇒ keine Nullstelle im Kreis)")
    print(f"   Mittel arg = {np.mean(ph - 2*np.pi*np.round(ph[M//2]/(2*np.pi))):+.3e}   Mittel ln| | = {np.mean(np.log(abs(G))):+.3e}")
    res[f"A_{name}"] = {"wind": wind, "mean_arg": float(np.mean(ph)), "mean_log": float(np.mean(np.log(abs(G))))}

# Φ auf einem Gitter, Extrapolation η → 0, zwei Datensätze
print("\n s     | Φ(s) (a1)   | Φ(s) (a2)   | |1+θ'| (a1)")
grid = [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98]
tabA = []
for s in grid:
    g1 = 1 + theta_d_extra(a1, s)[1]
    g2 = 1 + theta_d_extra(a2, s)[1]
    tabA.append([s, float(np.angle(g1)), float(np.angle(g2)), float(abs(g1))])
    print(f"{s:5.2f} | {np.angle(g1):+.8f} | {np.angle(g2):+.8f} | {abs(g1):.6f}")
res["Phi_grid"] = tabA

# ================= Teil B =================
print("\n== Teil B: reelle Achse, s + θ(s) = h(x+i0) ==")
tabB = []
for x in ["1.2", "1.5", "1.8", "2.0", "2.3", "2.6"]:
    hx, hpx = hp_real(mp.mpf(x))
    hx, hpx = complex(hx), complex(hpx)
    best = None
    for j in range(-3, 4):
        target = hx + j * complex(per)
        s = target - theta_d_extra(a1, 0.5)[0]
        s = complex(s.real % 1, s.imag) if False else s
        for _ in range(30):
            if not (0.01 < s.real % 1 < 0.99):
                break
            th, thp = theta_d_extra(a1, s.real)
            # Newton nur in reeller Richtung zulassen: komplexes s über Taylor
            f = s + th + 1j * s.imag * thp - target
            s = s - f / (1 + thp)
        if 0.01 < s.real % 1 < 0.99 and (best is None or abs(s.imag) < abs(best[0].imag)):
            best = (s, j)
    s, j = best
    th, thp = theta_d_extra(a1, s.real)
    sexpp = (1 + thp) / hpx
    tabB.append([x, s.real, s.imag, j, sexpp.real, sexpp.imag, float(np.angle(hpx)), float(np.angle(1 + thp))])
    print(f"x={x}: s={s.real:.10f} Im s={s.imag:+.2e} (Zweig j={j}); sexp'(s)=(1+θ')/h' = {sexpp.real:.8f} {sexpp.imag:+.2e}i;"
          f"  arg h'(x+i0)={np.angle(hpx):+.7f}  arg(1+θ'(s))={np.angle(1+thp):+.7f}")
res["B"] = tabB

# ================= Teil C =================
print("\n== Teil C: a_n gegen h'(n+i0) und −h'(c−ln n)/n ==")
tabC = []
for n in [2, 4, 16, 64, 256, 1024, 2048, 4096]:
    _, hpn = hp_real(mp.mpf(n))
    _, hpv = h_hp((c - mp.log(n)) / lb)
    lead = -hpv / lb / n
    an = a1[n]
    tabC.append([n, float(np.angle(an)), float(mp.arg(lead)), float(mp.arg(hpn)), float(abs(an / complex(hpn)))])
    print(f"n={n:5d}  arg a_n={np.angle(an):+.6f}  arg führend={float(mp.arg(lead)):+.6f}  arg h'(n)={float(mp.arg(hpn)):+.6f}"
          f"  |a_n/h'(n)|={abs(an/complex(hpn)):.5f}")
res["C"] = tabC
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / f"p1_{B}.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
