"""k59b — Kontrolle von k59 in hoher Genauigkeit: Wachstumsrate ‖L^k v‖^{1/k} der linearisierten Dyson-Abbildung, Variante „c-rü“
(nur Kompositionskanal Ẽ_c, Eichung: letzte Zeile von X gegeben). arb-Arithmetik, Potenzmethode mit Normierung.
Aufruf: python k59b_potenz.py n schritte   (z. B. 40 400)
"""
import math
import os
import random
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, modell  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def main():
    n = int(sys.argv[1]); K = int(sys.argv[2])
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    X = Jr.inv()
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    S = arb_mat(n + 1, n + 1)
    for k in range(1, n + 1):
        S[k, k - 1] = rho
    E = Jr * S - S * Jr
    fk, _ = modell("kubisch", n)
    f = arb_poly(fk[:4])
    pw = arb_poly([1])
    for _ in range(n + 1):
        pw = arb_poly([pw[k] if k <= pw.degree() else arb(0) for k in range(n + 1)]) * f
    for k in range(n + 1):
        ckk = pw[k] if k <= pw.degree() else arb(0)
        E[k, n] = E[k, n] + ckk * rho ** k / rho ** n
    A = E * X
    B = X * E
    random.seed(1)
    V = arb_mat(n + 1, n + 1)
    for i in range(n):
        for j in range(n + 1):
            V[i, j] = arb(random.uniform(-1, 1))
    logs = [0.0]
    for it in range(1, K + 1):
        R = V * A + B * V
        Y = arb_mat(n + 1, n + 1)
        for i in range(n, 0, -1):
            for j in range(n):
                Y[i - 1, j] = Y[i, j + 1] + R[i, j] / rho
            Y[i - 1, n] = R[i, n] / rho
        nrm = math.sqrt(sum(float(Y[i, j].mid()) ** 2 for i in range(n + 1) for j in range(n + 1)))
        logs.append(logs[-1] + math.log(nrm))
        V = arb_mat(n + 1, n + 1)
        for i in range(n + 1):
            for j in range(n + 1):
                V[i, j] = arb(Y[i, j].mid()) / nrm
        if it % 50 == 0:
            rate = math.exp((logs[it] - logs[it - 50]) / 50)
            print("n = %d, k = %d: (‖L^k v‖/‖L^{k−50} v‖)^{1/50} = %.4f  (%.0f s)" % (n, it, rate, time.time() - t0))
    print("n = %d: mittlere Rate k = 100 … %d: %.4f" % (n, K, math.exp((logs[K] - logs[100]) / (K - 100))))
    # Lage des dominanten Vektors: Masse nach Indexbereichen (i = Lösungsindex, j = Datenindex) und Kern-Werte δK(z,t) = ε_zᵀ V ẽ_t
    import cmath
    Vm = [[float(V[i, j].mid()) for j in range(n + 1)] for i in range(n + 1)]
    tot = sum(v * v for row in Vm for v in row)
    bins = [0, n // 4, n // 2, 3 * n // 4, n + 1]
    print("  Masse |V_ij|² nach (i-Viertel, j-Viertel):")
    for a_ in range(4):
        print("   i ∈ [%3d,%3d): " % (bins[a_], bins[a_ + 1]) + "  ".join("%.3f" % (sum(Vm[i][j] ** 2 for i in range(bins[a_], bins[a_ + 1]) for j in range(bins[b_], bins[b_ + 1])) / tot) for b_ in range(4)))
    Lc = complex(float(arb(L_RE).mid()), float(arb(L_IM).mid())); r = abs(Lc)
    zs = {"0": 0j, "0.5ρ": 0.5 * r, "L−0.02": Lc * (1 - 0.02 / r), "L−0.2": Lc * (1 - 0.2 / r), "Steg −ρ": -r, "Γout ρ": r}
    tt = {"1.02L": 1.02 * Lc, "L+0.2": Lc * (1 + 0.2 / r), "−1.02ρ": -1.02 * r, "1.02ρ": 1.02 * r, "2ρ": 2 * r}
    def dK(z, t_):
        s_ = 0j
        for i in range(n + 1):
            zi = (z / r) ** i
            for j in range(n + 1):
                s_ += zi * Vm[i][j] * (r / t_) ** j
        return abs(s_)
    print("  |δK(z,t)| (Zeilen z, Spalten t: " + ", ".join(tt) + ")")
    for zn, z in zs.items():
        print("   %-8s " % zn + "  ".join("%.2e" % dK(z, t_) for t_ in tt.values()))


if __name__ == "__main__":
    main()
