"""Strang A, Runde 2: Anpassungen an die Daten von Π_N f_K (N = 100 … 1000).

(a) c_N: freier komplexer Exponent mit/ohne 1/N-Glied, Schwingungsglieder (Frequenzen arg L und 2 arg L);
(b) K_N = c_N (δ_{−1})_N/(N u_N): Mittel über Fünfergruppen, Extrapolation K_∞ mit 1/N bzw. N^{−0.117};
(c) |P_N(L) − L| und max_Rand|P_N − f_K|: lokale Exponenten und Fit;
(d) Lage des Randmaximums.
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

HIER = Path(__file__).resolve().parent
DATEN = HIER / "daten"
mp.mp.dps = 30
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
th = mp.arg(L)
A = json.loads((DATEN / "lin_auswertung.json").read_text(encoding="utf-8"))
R = json.loads((DATEN / "lin_kreisfehler.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
Ns = sorted(int(k) for k in A)
c = {N: cx(A[str(N)]["c"]) for N in Ns}
K = {N: cx(A[str(N)]["K"]) for N in Ns}
PL = {N: cx(A[str(N)]["PL"]) for N in Ns}


def unwrap(vals):
    out = []
    for v in vals:
        if out:
            while v - out[-1] > mp.pi:
                v -= 2 * mp.pi
            while v - out[-1] < -mp.pi:
                v += 2 * mp.pi
        out.append(v)
    return out


def fit(Nsel, werte, basis):
    """komplexe Grösse w_N: log|w| und arg w (stetig) getrennt nach basis(N) anpassen."""
    Am = mp.matrix([basis(N) for N in Nsel])
    br = mp.qr_solve(Am, mp.matrix([mp.log(abs(werte[N])) for N in Nsel]))[0]
    bi = mp.qr_solve(Am, mp.matrix(unwrap([mp.arg(werte[N]) for N in Nsel])))[0]
    return br, bi


osz = lambda N: [mp.cos(N * th) / N, mp.sin(N * th) / N, mp.cos(2 * N * th) / N, mp.sin(2 * N * th) / N]  # noqa: E731
print("(a) c_N ~ N^β: freie Anpassung (Soll 1 − 2γ = " + mp.nstr(1 - 2 * G, 6) + ")")
for lo, hi in ((100, 400), (300, 1000), (500, 1000), (100, 1000)):
    sel = [N for N in Ns if lo <= N <= hi]
    for name, bas in (("ohne 1/N", lambda N: [1, mp.log(N)] + osz(N)), ("mit 1/N", lambda N: [1, mp.log(N), mp.mpf(1) / N] + osz(N)),
                      ("mit 1/N, 1/N²", lambda N: [1, mp.log(N), mp.mpf(1) / N, mp.mpf(1) / N ** 2] + osz(N))):
        br, bi = fit(sel, c, bas)
        print(f"   N ∈ [{lo},{hi}] ({len(sel)} Werte), {name:>14}: β = {mp.nstr(br[1], 6)} + {mp.nstr(bi[1], 6)} i")
print("\n(b) K_N: Fünfergruppen (Mittel über N0−4 … N0) und Einzelwerte")
gruppen = {}
for N0 in (300, 500, 700, 1000):
    g = [K[N] for N in range(N0 - 4, N0 + 1) if N in K]
    if len(g) == 5:
        gruppen[N0] = mp.fsum(g) / 5
        print(f"   N0 = {N0}: K = {mp.nstr(gruppen[N0], 7)}  (Streuung Re {mp.nstr(max(mp.re(x) for x in g) - min(mp.re(x) for x in g), 2)}, "
              f"Im {mp.nstr(max(mp.im(x) for x in g) - min(mp.im(x) for x in g), 2)})")
# Anpassung K_N = K∞ + b/N + Schwingung über alle N ≥ 200
sel = [N for N in Ns if N >= 200]
for name, fn in (("1/N", lambda N: mp.mpf(1) / N), ("N^-0.117", lambda N: mp.power(N, -0.117)), ("N^-0.5", lambda N: mp.power(N, -0.5))):
    Am = mp.matrix([[1, fn(N)] + [x * N for x in osz(N)] for N in sel])
    xr = mp.qr_solve(Am, mp.matrix([mp.re(K[N]) for N in sel]))
    xi = mp.qr_solve(Am, mp.matrix([mp.im(K[N]) for N in sel]))
    print(f"   K_N = K∞ + b·{name} + Schwingung: K∞ = {mp.nstr(xr[0][0] + 1j * xi[0][0], 6)}, b = {mp.nstr(xr[0][1] + 1j * xi[0][1], 4)}, "
          f"Rest {mp.nstr(xr[1], 2)}/{mp.nstr(xi[1], 2)}")
print("\n(c) |P_N(L) − L| und max_Rand|P_N − f_K|: lokale Exponenten")
prev = None
for N in Ns:
    if N % 100 != 0 and N not in (150, 250):
        continue
    a1, a2 = abs(PL[N]), R[str(N)]["maxE"]
    if prev:
        N1, b1, b2 = prev
        print(f"   {N1}→{N}: |P(L)−L| {mp.nstr(-mp.log(a1 / b1) / mp.log(mp.mpf(N) / N1), 4)}, max_Rand {mp.nstr(-mp.log(a2 / b2) / mp.log(mp.mpf(N) / N1), 4)};"
              f"  Werte {mp.nstr(a1, 4)}, {a2:.4g}; Ort des Randmaximums (θ − arg L)·N = {R[str(N)]['ort_N']:.3g}")
    prev = (N, a1, a2)
sel = [N for N in Ns if N >= 300]
br, bi = fit(sel, PL, lambda N: [1, mp.log(N), mp.mpf(1) / N] + osz(N))
print(f"   Fit N ≥ 300: P_N(L) − L ~ N^({mp.nstr(br[1], 5)} + {mp.nstr(bi[1], 5)} i)  [−γ = {mp.nstr(-G, 5)}; 5.447/1.1325 = 4.810]")
br, bi = fit(sel, PL, lambda N: [1, mp.log(N)] + osz(N))
print(f"   Fit N ≥ 300 ohne 1/N: N^({mp.nstr(br[1], 5)} + {mp.nstr(bi[1], 5)} i)")
print("\n(d) max|P_N| auf D (global, nahe L):", ", ".join(f"{N}: {R[str(N)]['maxP_global']:.6f}/{R[str(N)]['maxP_nahL']:.6f}"
                                                    for N in (100, 300, 500, 700, 1000) if str(N) in R))
