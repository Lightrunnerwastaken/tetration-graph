"""s1: Knesers Erzeuger J_K = sexp'∘slog nahe am Fixpunkt L (Basis e).

J_K(w) = J_L(w) * (1 + θ'(s)),  s = slog_K(w) löst s + θ(s) = h(w),  J_L = 1/h' (regulärer Erzeuger).
Prüft Satz T4: J_K - J_L ≈ 2πi a_1 e^{2πi s} J_L  ~  C (w-L)^{1+2πi/L}  (Exponent-Realteil 1 + 2π Im L/|L|^2).
Gegenprobe: zwei Genauigkeiten (dps 40 / 60, K 50 / 70) und Motorwert sexp an einem Punkt.
Aufruf: python s1_erzeuger_bei_L.py
"""
import sys
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib10 as lib


def slog_K(T, B, w, s0=None):
    hw, dh = B.h_mit_ableitung(w)
    P = 2j * mp.pi / B.L          # Periode von h (h ist nur modulo 2πi/L bestimmt)
    s = hw - T.a[0]
    k = mp.floor((mp.re(s) + 1) / mp.re(P))
    s = s - k * P                 # Zweig mit Re s in [-1, -1+4.448): fest gewählter Kneser-Zweig
    hw = s + T.a[0]
    for _ in range(60):
        th, th1 = T.theta(s, 1)
        ds = (s + th - hw) / (1 + th1)
        s -= ds
        if abs(ds) < mp.mpf(10) ** (-mp.mp.dps + 3):
            break
    th, th1 = T.theta(s, 1)
    assert abs(s + th - hw) < mp.mpf(10) ** (-mp.mp.dps + 8), "Newton"
    JL = 1 / dh
    return s, JL, JL * (1 + th1)


def lauf(dps, K):
    mp.mp.dps = dps
    B = lib.Basis(K=K)
    a = lib.lade_theta()
    T = lib.ThetaReihe(a, B)
    L = B.L
    alpha = 1 + 2j * mp.pi / L
    out = []
    for phi in [mp.mpf(0), mp.pi / 3, 2 * mp.pi / 3, -mp.pi / 2]:
        for r in [mp.mpf("1e-1"), mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4")]:
            w = L + r * mp.expj(phi)
            s, JL, JK = slog_K(T, B, w)
            D = JK - JL
            # Vorhersage aus dem ersten θ-Glied
            pred = JL * 2j * mp.pi * a[1] * mp.expjpi(2 * s)
            out.append((phi, r, D, pred, s))
    return L, alpha, a, out


if __name__ == "__main__":
    L, alpha, a, o1 = lauf(40, 50)
    _, _, _, o2 = lauf(60, 70)
    mp.mp.dps = 30
    print("L =", mp.nstr(L, 20), " Exponent 1+2πi/L =", mp.nstr(alpha, 15))
    print("a_0 =", mp.nstr(a[0], 15), " a_1 =", mp.nstr(a[1], 15))
    print("phi      r       |D|=|J_K-J_L|     |D/pred - 1|    Abweichung dps40/60")
    for (phi, r, D, pred, s), (_, _, D2, _, _) in zip(o1, o2):
        print(f"{float(phi):6.3f} {float(r):8.0e}  {mp.nstr(abs(D), 8):>14}  {mp.nstr(abs(D/pred-1), 4):>12}  {mp.nstr(abs(D-D2)/abs(D), 3):>10}")
    # lokaler Exponent: log|D(r)|/log r Steigung
    print("\nlokale Steigung d log|D| / d log r (Vorhersage Re α =", mp.nstr(mp.re(alpha), 12), ")")
    for i in range(0, len(o1), 4):
        rows = o1[i:i + 4]
        for j in range(3):
            sl = (mp.log(abs(rows[j + 1][2])) - mp.log(abs(rows[j][2]))) / (mp.log(rows[j + 1][1]) - mp.log(rows[j][1]))
            print(f"  phi={float(rows[j][0]):6.3f}  r {float(rows[j][1]):.0e}->{float(rows[j+1][1]):.0e}: {mp.nstr(sl, 12)}")
    # Koeffizient: D ≈ C (w-L)^α ; C = 2πi a_1 e^{-2πi a_0} L * (Normierung σ'(L)=1)
    C_pred = 2j * mp.pi * a[1] * mp.exp(-2j * mp.pi * a[0]) * L
    print("\nC_vorhergesagt = 2πi a_1 e^{-2πi a_0} L =", mp.nstr(C_pred, 15))
    for (phi, r, D, pred, s) in o1:
        if r == mp.mpf("1e-4"):
            w_L = r * mp.expj(phi)
            print(f"  phi={float(phi):6.3f}: D/(w-L)^α =", mp.nstr(D / mp.exp(alpha * mp.log(w_L)), 15))
