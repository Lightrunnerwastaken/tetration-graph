"""M28 — Kante (a″)_{αN} im radialen Modell: Wachstum des Überlaufblocks Q_{(1+α)N}J_∞P_N und σ_min der Rechteck-Abschnitte.
Radiales Modell: f quadratisch ⇒ J_∞P_N ⊂ P_{2N} (alle Zeilen exakt vorhanden). arb (prec gross), Spaltennormen und σ_min über (RᵀR)⁻¹.
Aufruf: python m28_kante_ueberlauf.py "N1,N2" "a1,a2,..."
"""
import sys

from flint import arb, arb_mat, ctx

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "100,200").split(",")]
ALS = [float(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "0.01,0.03,0.05,0.1,0.3,1").split(",")]
ctx.prec = 800
eps = arb("0.31382")
th = arb("76.618") * arb.pi() / 180
c0 = 1 - 2 * eps * th.cos()
fco = [eps, c0, eps]                  # f = ε + c0 z + ε z²
# m = f'∘f = c0 + 2ε f
mco = [c0 + 2 * eps * eps, 2 * eps * c0, 2 * eps * eps]


def mul(a, b):
    out = [arb(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


for N in NS:
    R = 2 * N
    pot = [[arb(1)]]
    for j in range(1, N + 1):
        pot.append(mul(pot[-1], fco))
    J = arb_mat(R + 1, N + 1)
    for j in range(N + 1):
        for k, c in enumerate(pot[j]):
            J[k, j] += c
        for i, c in enumerate(mco):
            if j + i <= R:
                J[j + i, j] += c
    for al in ALS:
        M = min(R, N + max(0, int(round(al * N))))
        # Überlaufblock: Zeilen M+1..R; grösste Spaltennorm und Frobenius
        maxcol = arb(0)
        fro = arb(0)
        for j in range(N + 1):
            s = arb(0)
            for k in range(M + 1, R + 1):
                s += J[k, j] ** 2
            fro += s
            if s.upper() > maxcol.upper():
                maxcol = s
        # σ_min des Rechteck-Abschnitts T_M J
        T = arb_mat(M + 1, N + 1)
        for k in range(M + 1):
            for j in range(N + 1):
                T[k, j] = J[k, j]
        G = T.transpose() * T
        try:
            Gi = G.inv()
            # λ_max(G⁻¹) per Potenziteration (mid), σ_min = λ_max^{-1/2}
            import mpmath as mp
            mp.mp.dps = 60
            A = mp.matrix([[mp.mpf(Gi[i, j].mid().str(60, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
            ev = mp.eigsy(A, eigvals_only=True)
            smin = 1 / mp.sqrt(max(ev))
        except Exception as e:
            smin = f"nan ({e})"
        print(f"N = {N}, α = {al}: Zeilen bis {M}: max Spaltennorm Überlauf = {float(maxcol.sqrt().mid()):.3e}, "
              f"Frobenius = {float(fro.sqrt().mid()):.3e};  σ_min(T_{{{M}}}J_∞|P_N) = {float(smin) if not isinstance(smin, str) else smin}",
              flush=True)
