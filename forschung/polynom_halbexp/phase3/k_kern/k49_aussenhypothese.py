"""k49 — Vorprüfung der Aussenhypothese (AH) in Kernform (BERICHT §28):
  Q_n = max_{z ∈ ∂K_{r0}} max_{t ∈ Γ*} |ΔK_n(z,t)|,   ΔK_n(z,t) = δ_n(z)Λ_n(t)/ℓ_nn  (Rang-1-Zuwachs des Eckkerns K_n(z,t) = ε_zᵀJ̃_n⁻¹ẽ_t)
  = |P_n| · R_P · R_A,   R_P = max_{∂K_{r0}}|δ_n/δ_n(L)|,   R_A = max_{Γ*}|Λ_n/Λ_n(L)|.
Γ* nahe L: zwei FESTE Log-Spiralen θ(r) = θ_j + 1.623°·ln(r/ρ) (Fenstermitte u_A = 60° bzw. −90°, arg Λ₀ = −71.46°, vgl. V §3.18), r ∈ [1e-7, r1];
fern (nur grob, V macht es richtig): |t| = ρ(1 + S) über Γin (Winkelabstand > 12° von L, L̄), |t| = ρ über Γout, S = 0.0166.
Zusätzlich: Maximum auf dem Kegelkreis |t − L| = A0·r_n (alle Winkel) — für das Maximumprinzip im Eckteil.
Aufruf: python k49_aussenhypothese.py 200,400,600 [u_A-Liste, Standard 60,-90] [r1, Standard 0.1]
"""
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    uAs = [float(s) for s in sys.argv[2].split(",")] if len(sys.argv) > 2 else [60.0, -90.0]
    r1 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        Jc = acb_mat(J)
        e = acb_mat(n + 1, 1); e[n, 0] = acb(1)
        c = Jc.solve(e, algorithm="approx")
        l = Jc.transpose().solve(e, algorithm="approx")
        mp.mp.dps = int(bits / 3.32) + 20
        L = mp.mpc(L_RE, L_IM); rho = abs(L); Lb = mp.conj(L)
        nd = mp.mp.dps

        def cv(x):
            return mp.mpc(mp.mpf(x.real.mid().str(nd, radius=False)), mp.mpf(x.imag.mid().str(nd, radius=False)))
        cm = [cv(c[k, 0]) for k in range(n + 1)]
        lm = [cv(l[k, 0]) for k in range(n + 1)]

        def delta(z):
            w = z / rho; s = mp.mpc(0)
            for k in range(n, -1, -1):
                s = s * w + cm[k]
            return s

        def Lam(t):
            w = rho / t; s = mp.mpc(0)
            for k in range(n, -1, -1):
                s = s * w + lm[k]
            return s
        dL = delta(L); LL = Lam(L); P = dL * LL / lm[n]
        rn = rho * mp.power(n, -0.8935896)
        # primal: Rand von K_{r0}
        out = []
        for r0 in [0.1, 0.05, 0.02]:
            best = (mp.mpf(0), None)
            M = 2400
            for j in range(M):
                ph = 2 * mp.pi * j / M
                z = rho * mp.expj(ph)
                if abs(z - L) < r0 or abs(z - Lb) < r0:
                    continue
                v = abs(delta(z) / dL)
                if v > best[0]:
                    best = (v, "Kreis %.1f°" % float(mp.degrees(ph)))
            for cen in (L, Lb):
                for j in range(720):
                    z = cen + r0 * mp.expj(2 * mp.pi * j / 720)
                    if abs(z) > rho:
                        continue
                    v = abs(delta(z) / dL)
                    if v > best[0]:
                        best = (v, "Bogen um %s, %.0f°" % ("L" if cen == L else "L̄", 360.0 * j / 720))
            out.append("r0=%.2f: R_P = %s (%s)" % (r0, mp.nstr(best[0], 4), best[1]))
        # adjungiert: Spiralen
        argL0 = -71.46
        spir = []
        Rspir = mp.mpf(0)
        for uA in uAs:
            th0 = uA - argL0
            mx = (mp.mpf(0), None)
            for i in range(400):
                lr = mp.log(1e-7) + (mp.log(r1) - mp.log(1e-7)) * i / 399
                r = mp.exp(lr)
                th = th0 + 1.623 * mp.log(r / rho)
                t = L + r * mp.expj(mp.radians(th))
                v = abs(Lam(t) / LL)
                if v > mx[0]:
                    mx = (v, "r = %s (A = %s), θ = %.1f°" % (mp.nstr(r, 3), mp.nstr(r / rn, 3), float(th)))
            spir.append("u_A=%+.0f°: max %s bei %s" % (uA, mp.nstr(mx[0], 4), mx[1]))
            Rspir = max(Rspir, mx[0])
        # grob fern
        S = 0.0166
        mxf = (mp.mpf(0), None)
        argLdeg = float(mp.degrees(mp.arg(L)))
        M = 2400
        for j in range(M):
            ph = 360.0 * j / M
            dist = min(abs(ph - argLdeg), abs(ph - (360 - argLdeg)))
            if dist < 12:
                continue
            gin = argLdeg < ph < 360 - argLdeg
            t = rho * (1 + (S if gin else 0)) * mp.expj(mp.radians(ph))
            v = abs(Lam(t) / LL)
            if v > mxf[0]:
                mxf = (v, "%.1f° (%s)" % (ph, "Γin" if gin else "Γout"))
        # Kegelkreis
        kk = []
        for A0 in (1, 2):
            mx = mp.mpf(0)
            for j in range(360):
                t = L + A0 * rn * mp.expj(2 * mp.pi * j / 360)
                mx = max(mx, abs(Lam(t) / LL))
            kk.append("A0=%d: %s" % (A0, mp.nstr(mx, 4)))
        print("n = %d (%.0f s): P_n = %s, |P_n| = %s, r_n = %s" % (n, time.time() - t0, mp.nstr(P, 6), mp.nstr(abs(P), 5), mp.nstr(rn, 3)), flush=True)
        print("   primal: " + " | ".join(out), flush=True)
        print("   adjungiert Spiralen: " + " | ".join(spir), flush=True)
        print("   adjungiert fern (grob): max %s bei %s" % (mp.nstr(mxf[0], 4), mxf[1]), flush=True)
        print("   adjungiert Kegelkreis max|Λ/Λ(L)|: " + " | ".join(kk), flush=True)


if __name__ == "__main__":
    main()
