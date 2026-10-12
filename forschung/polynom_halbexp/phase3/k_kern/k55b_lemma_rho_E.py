"""k55b — Vortest (b), richtige Fassung: Parametrix als Approximant auf E statt Taylor-Abschnitt (BERICHT §37). Kubisches Modell.
Zusammengesetzte Funktion c_n(x) = s_reg(x) + (C/d_n)·V((x − L)/d_n)  (äussere Lösung ohne Pol + inneres Profil statt Pol).
ũ_n := Kleinste-Quadrate-Polynom vom Grad n zu c_n auf Randpunkten von E (Vandermonde mit Arnoldi, double), Koeffizienten per FFT auf |x| = ρ(1 − 3/n).
E wie k54 mit (u1, u2) (Standard 80°, −110°: dann liegt f(Σ₂) ausserhalb des primalen Wachstumsfensters).
s_reg auf E: in D Taylor-Reihe (Galerkin J̃_{2n}⁻¹T g, g = J s_reg regulär); ausserhalb D über die Funktionalgleichung s_reg∘f + m·s_reg = g.
Dann Lemma ρ: ρ_{m,n} = (L⁻¹r̃)_m/(L⁻¹b)_m, r̃ = b − J̃_nũ (hohe Genauigkeit, LU ohne Pivot).
Aufruf: python k55b_lemma_rho_E.py n [u1 u2]
"""
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff, modell  # noqa: E402
import k54_geometrie_E as G  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def acb2mp(z, nd):
    return mp.mpc(mp.mpf(z.real.mid().str(nd, radius=False)), mp.mpf(z.imag.mid().str(nd, radius=False)))


def main():
    n = int(sys.argv[1])
    u1 = float(sys.argv[2]) if len(sys.argv) > 2 else 80.0
    u2 = float(sys.argv[3]) if len(sys.argv) > 3 else -110.0
    t0 = time.time()
    bits = int(1.3 * 2 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); Lb = mp.conj(L); rho = abs(L)
    al = L / (4 * mp.pi); a = 1 + 1j * al; mu = mp.sqrt(L)
    cs = -1 / (mu + 1 / mu)
    bP = 0.5 + 1j * al
    lamP = mp.mpc("0.603424", "1.165020"); lamPt = -lamP * rho ** 2 / L ** 2
    Ck = -lamPt * mp.gamma(0.5 + 1j * al) / ((mu + 1 / mu) * mp.gamma(0.5 + 2j * al))
    dn = rho * mp.power(n, -a)
    N = 2 * n
    # ---- g und s_reg (wie k55)
    ck = [mp.mpf(x.mid().str(nd, radius=False)) for x in kubik_koeff()]
    alpha = ck[2] + 2 * mp.re(L) * ck[3]; beta = ck[3]
    Q = [1 - Lb * alpha, alpha - Lb * beta, beta]
    invQ = [mp.mpc(0)] * (N + 1); invQ[0] = 1 / Q[0]
    for k in range(1, N + 1):
        invQ[k] = -(Q[1] * invQ[k - 1] + (Q[2] * invQ[k - 2] if k >= 2 else 0)) / Q[0]
    invxL = [-1 / L ** (k + 1) for k in range(N + 1)]
    B = [mp.fsum(invxL[j] * invQ[k - j] for j in range(k + 1)) for k in range(N + 1)]
    fk, mk = modell("kubisch", 12)
    m = [mp.mpf(x.mid().str(nd, radius=False)) for x in mk[:7]]
    mL = mp.polyval(m[::-1], L)
    quot = []; acc = mp.mpc(0)
    for c in m[::-1][:-1]:
        acc = acc * L + c; quot.append(acc)
    quot = quot[::-1]
    Cm = [(quot[k] if k < len(quot) else 0) + mL * invxL[k] for k in range(N + 1)]
    g = [1 / L ** (k + 1) - cs * (B[k] + Cm[k]) for k in range(N + 1)]
    J2 = acb_mat(J_tilde("kubisch", N, bits=bits)); ctx.prec = bits
    gv = acb_mat(N + 1, 1)
    for k in range(N + 1):
        v = g[k] * rho ** k
        gv[k, 0] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    sreg = J2.solve(gv, algorithm="approx")
    sr_full = [complex(acb2mp(sreg[k, 0], 30)) / float(rho) ** k for k in range(N + 1)]   # echte Taylor-Koeffizienten (double)
    print("n = %d, Fenster (%.0f°, %.0f°): s_reg fertig (%.0f s); Gewichtung %s" % (n, u1, u2, time.time() - t0, os.environ.get("K55_GEWICHT", "0")))
    # ---- Hilfsfunktionen in double
    Lc = complex(L); rhoc = float(rho); csc = complex(cs)
    fc = np.array([complex(x) for x in ck])
    mcoef = np.array([float(x) for x in m])
    def f(x): return fc[0] + x * (fc[1] + x * (fc[2] + x * fc[3]))
    def mfun(x): return np.polyval(mcoef[::-1], x)
    def gfun(x): return 1 / (Lc - x) - csc * (1 / (f(x) - Lc) + mfun(x) / (x - Lc))
    srp = np.array(sr_full[::-1])
    def s_in(x): return np.polyval(srp, x)                 # Taylor in D (und auf dem Kreis)
    def s_reg(x):
        x = np.asarray(x, dtype=complex)
        out = np.empty_like(x)
        inside = np.abs(x) <= rhoc * (1 + 1e-9)
        out[inside] = s_in(x[inside])
        xo = x[~inside]
        fxo = f(xo)
        if len(xo) and np.max(np.abs(fxo)) >= rhoc:
            raise RuntimeError("Funktionalgleichung: f(x) nicht in D (max |f|/ρ = %.4f)" % (np.max(np.abs(fxo)) / rhoc))
        out[~inside] = (gfun(xo) - s_in(fxo)) / mfun(xo)   # f(x) ∈ D (Linse, Σ₁)
        return out
    # ---- V (mpmath, Koeffizienten vorab)
    umax = 1.32 * 4.0 / abs(complex(dn))
    Dmax = 60 + int(0.45 * umax ** 1.12)
    Kmax = int(3.5 * umax ** 1.12) + 200
    mp.mp.dps = Dmax
    gb = mp.gamma(bP)
    cV = []
    for k in range(Kmax):
        cV.append(mp.gamma(bP - 1j * al * k) / (mp.factorial(k) * gb) / (k + 1))
    print("   V-Koeffizienten: %d Terme bei %d Stellen (|u|max ≈ %.0f)" % (Kmax, Dmax, umax))
    def V(z):
        w = -lamPt * mp.mpc(z)
        aw = abs(w)
        dps = min(Dmax, 40 + int(0.45 * float(aw) ** 1.12))
        with mp.workdps(dps):
            s_ = mp.mpc(0); p = mp.mpc(1); big = mp.mpf(0)
            for k in range(Kmax):
                t = cV[k] * p
                s_ += t; at = abs(t); big = max(big, at)
                if k > 20 and at < mp.mpf(10) ** (-25) * (abs(s_) + 1e-300) and at < big * mp.mpf(10) ** (-dps + 25):
                    break
                p *= w
            return complex(s_)
    dnc = complex(dn); Ckc = complex(Ck)
    def cfun(x):
        return s_reg(np.array([x]))[0] + Ckc / dnc * V((x - Lc) / dnc)
    # ---- Randpunkte von E (aus k54): Γin-Seite von 𝒞 direkt, Γout-Seite als f-Bild
    parts, names = G.kontur(u1, u2, 0.03, 0.05, m=300)
    rmin = 0.1 * abs(complex(dn))
    def keep(z):
        z = np.asarray(z)
        return z[np.abs(z - Lc) > rmin][::2] if len(z) else z
    gin_side = [parts[0], parts[1], parts[2], parts[3], parts[4]]          # Σ1, V1, Linse, V1b, Σ1b
    gout_side = [parts[5], parts[6], parts[7], parts[8], parts[9]]         # Σ2b, V2b, Γout, V2, Σ2
    gew = os.environ.get("K55_GEWICHT", "0") == "1"
    T_all = np.concatenate([keep(p) for p in gin_side + gout_side])
    T_all = T_all[(np.abs(T_all - Lc) > rmin) & (np.abs(T_all - Lc.conjugate()) > rmin)]
    nahe = (np.abs(T_all - Lc) <= 0.05) | (np.abs(T_all - Lc.conjugate()) <= 0.05)
    aussenD = np.abs(T_all) > rhoc * (1 + 1e-9)
    wA = np.ones(len(T_all))
    if gew:
        wA[~nahe & aussenD] = math.exp(-0.013 * n)       # Linse
        wA[~nahe & ~aussenD] = math.exp(-0.037 * n)      # Γout-Bogen
    fT = f(T_all)
    srT = s_reg(T_all)                                   # in D̄ Taylor, ausserhalb über Funktionalgleichung (f(t) ∈ D)
    # s_reg(f(t)): falls f(t) ∈ D Taylor, sonst g(t) − m(t)s_reg(t)
    in_fT = np.abs(fT) <= rhoc * (1 + 1e-9)
    srfT = np.where(in_fT, s_in(np.where(in_fT, fT, 0)), gfun(T_all) - mfun(T_all) * srT)
    X = np.concatenate([T_all, fT]); SR = np.concatenate([srT, srfT]); Wt = np.concatenate([wA, wA])
    xin = T_all
    t1 = time.time()
    Vv = np.array([V((x - Lc) / dnc) for x in X])
    Cvals = SR + Ckc / dnc * Vv
    print("   %d Randpunkte, V ausgewertet (%.0f s); max|c_n| auf ∂E = %.3g, max|V-Teil| = %.3g" % (len(X), time.time() - t1, np.max(np.abs(Cvals)), np.max(np.abs(Ckc / dnc * Vv))))
    # ---- LSQ-Fit (Vandermonde mit Arnoldi)
    zc = 0.0; sc = 1.6 * rhoc
    xs = (X - zc) / sc
    Mp = len(X)
    Qm = np.zeros((Mp, n + 1), dtype=complex); H = np.zeros((n + 1, n), dtype=complex); Qm[:, 0] = 1
    for k in range(n):
        v = xs * Qm[:, k]
        for j in range(k + 1):
            H[j, k] = np.vdot(Qm[:, j], v) / Mp; v = v - H[j, k] * Qm[:, j]
        H[k + 1, k] = np.linalg.norm(v) / math.sqrt(Mp); Qm[:, k + 1] = v / H[k + 1, k]
    coef = np.linalg.lstsq(Qm * Wt[:, None], Cvals * Wt, rcond=None)[0]
    def qeval(z):
        z = (np.asarray(z) - zc) / sc
        W = np.zeros((len(z), n + 1), dtype=complex); W[:, 0] = 1
        for k in range(n):
            v = z * W[:, k]
            for j in range(k + 1):
                v = v - H[j, k] * W[:, j]
            W[:, k + 1] = v / H[k + 1, k]
        return W @ coef
    fiterr = np.max(np.abs(qeval(X) - Cvals))
    # Diagnose je Randstück
    err = np.abs(qeval(X) - Cvals)
    nin = len(xin)
    segs = [("Konturpunkte t", slice(0, nin)), ("Bilder f(t)", slice(nin, len(X)))]
    print("      gewichtet: %s; gewichteter max-Fehler = %.3g" % (gew, np.max(err * Wt)))
    for nm, sl in segs:
        e_ = err[sl]; c_ = np.abs(Cvals[sl]); v_ = np.abs(Ckc / dnc * Vv[sl]); s_ = np.abs(SR[sl])
        i = np.argmax(e_)
        print("      %s: max Fit-Fehler %.3g bei x = %s (|x−L| = %.3g); max|c| = %.3g, max|V-Teil| = %.3g (bei x = %s), max|s_reg| = %.3g"
              % (nm, e_[i], np.round(X[sl][i], 4), abs(X[sl][i] - Lc), c_.max(), v_.max(), np.round(X[sl][np.argmax(v_)], 4), s_.max()))
    # Fehler im Kegelbereich und in D (Kontrolle)
    # Monomkoeffizienten per FFT auf |x| = Rc
    Rc = rhoc * (1 - 3.0 / n)
    Mf = 8 * n
    xx = Rc * np.exp(2j * np.pi * np.arange(Mf) / Mf)
    qv = qeval(xx)
    ak = np.fft.fft(qv) / Mf
    utc = [ak[k] / Rc ** k * rhoc ** k for k in range(n + 1)]          # ρ-skaliert
    print("   Fit: max|q − c_n| auf ∂E = %.3g (relativ %.3g); FFT-Rest |a_k| für k > n: %.3g" % (fiterr, fiterr / np.max(np.abs(Cvals)), np.max(np.abs(ak[n + 1:Mf // 2]))))
    # ---- Lemma ρ in hoher Genauigkeit
    mp.mp.dps = nd
    Jn = acb_mat(J_tilde("kubisch", n, bits=bits)); ctx.prec = bits
    Jm = [[acb2mp(Jn[i, j], nd) for j in range(n + 1)] for i in range(n + 1)]
    bvec = [rho ** k / L ** (k + 1) for k in range(n + 1)]
    bv = acb_mat(n + 1, 1)
    for k in range(n + 1):
        bv[k, 0] = acb(arb(mp.nstr(mp.re(bvec[k]), nd)), arb(mp.nstr(mp.im(bvec[k]), nd)))
    ue = Jn.solve(bv, algorithm="approx")
    uex = [acb2mp(ue[k, 0], nd) for k in range(n + 1)]
    Aw = [row[:] for row in Jm]
    Lf = [[mp.mpc(0)] * (n + 1) for _ in range(n + 1)]
    for k in range(n + 1):
        Lf[k][k] = mp.mpc(1); piv = Aw[k][k]; rowk = Aw[k]
        for i in range(k + 1, n + 1):
            lik = Aw[i][k] / piv; Lf[i][k] = lik; rowi = Aw[i]
            for j in range(k + 1, n + 1):
                rowi[j] -= lik * rowk[j]
    def fwd(yv):
        z = [mp.mpc(0)] * (n + 1)
        for i in range(n + 1):
            s_ = yv[i]; Li = Lf[i]
            for j in range(i):
                s_ -= Li[j] * z[j]
            z[i] = s_
        return z
    yb = fwd(bvec)
    epsL = [(L / rho) ** k for k in range(n + 1)]
    uL = mp.fsum(epsL[k] * uex[k] for k in range(n + 1))
    def auswerten(name, ut):
        r = [bvec[i] - mp.fsum(Jm[i][j] * ut[j] for j in range(n + 1)) for i in range(n + 1)]
        yr = fwd(r)
        ab = [abs(yr[i] / yb[i]) for i in range(n + 1)]
        mx = max(range(n + 1), key=lambda i: ab[i])
        prof = " ".join("%.2f:%s" % (q, mp.nstr(ab[min(n, int(round(q * n)))], 3)) for q in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0))
        eL = mp.fsum(epsL[k] * (uex[k] - ut[k]) for k in range(n + 1))
        print("   %s: max_m|ρ_{m,n}| = %s bei m/n = %.3f | Profil: %s | |e_n|/|u_n(L)| = %s" % (name, mp.nstr(ab[mx], 4), mx / n, prof, mp.nstr(abs(eL) / abs(uL), 4)), flush=True)
    auswerten("Kontrolle exakt", uex)
    ut = [mp.mpc(complex(u)) for u in utc]
    auswerten("P3 (Approximant auf E)", ut)
    print("fertig (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
