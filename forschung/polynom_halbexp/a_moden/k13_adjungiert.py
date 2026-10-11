"""Strang A: adjungierte Mode. c_N^lin = ⟨ℓ, d⟩, d = −J⁻¹F(q_N) ⇒ c_N^lin = −⟨m_N, F(q_N)⟩, m_N = J^{−T} ℓ.

ℓ = Mode-Funktional (Koeffizient von e^{−2πiw} in θ = (p − f)(T(w))/T'(w+½), Mittel über w ∈ [−1.6, −0.6), y = 0, 48 Punkte).
Frage: welche Gleichungen k (Komponenten des Konsistenzfehlers F(q_N)) speisen die gefährliche Mode, wie stark löschen sie sich aus?
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, arb, arb_mat, arb_poly, ctx  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 60
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
qs = json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]
q = [mp.mpf(x) for x in qs]
am = json.loads((DATEN / "amplitude_cN.json").read_text(encoding="utf-8"))
kd = json.loads((DATEN / "kreis_r1.34_M4096.json").read_text(encoding="utf-8"))
M = kd["M"]
Ah = [mp.mpc(*x) for x in kd["A"]]
Av = Ah + [mp.conj(x) for x in Ah[1:-1][::-1]]
ctx.prec = 260
dd = acb.dft([acb(arb(mp.nstr(mp.re(x), 58)), arb(mp.nstr(mp.im(x), 58))) for x in Av])
Ak = [mp.mpf(dd[k].real.mid().str(58, radius=False)) / M / mp.mpf("1.34") ** k for k in range(1200)]


def reihe(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def abl(c, z):
    return reihe([k * c[k] for k in range(1, len(c))], z)


m = 48
ws = [mp.mpf("-1.6") + mp.mpf(j) / m for j in range(m)]
gp = FatouGP(dps=57, fatou_gp="fork")
zs = [mp.re(z) for z in gp.session("exp(1)").sexp_batch(ws)]
gp.close()
Tp = [mp.re(abl(q, z) / abl(Ak, z)) for z in zs]


def FJ(a, N):
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    g, fpf = arb_poly([0]), arb_poly([0])
    for mm in range(N + 1):
        g += a[mm] * P[mm]
        if mm:
            fpf += mm * a[mm] * P[mm - 1]
    F, J, fak = arb_mat(N + 1, 1), arb_mat(N + 1, N + 1), arb(1)
    for k in range(N + 1):
        fak = fak * k if k else fak
        F[k, 0] = g[k] - 1 / fak
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


t0 = time.time()
for N in (60, 100, 150):
    ctx.prec = 300 + 6 * N
    a = [arb(x) for x in qs[:N + 1]]
    F, J = FJ(a, N)
    # ℓ_k = (1/m) Σ_j z_j^k e^{2πi w_j} / T'_j  (komplex) → Real- und Imaginärteil getrennt
    lr = arb_mat(N + 1, 1)
    li = arb_mat(N + 1, 1)
    for k in range(N + 1):
        s = mp.fsum(zs[j] ** k * mp.expjpi(2 * ws[j]) / Tp[j] for j in range(m)) / m
        lr[k, 0] = arb(mp.nstr(mp.re(s), 55))
        li[k, 0] = arb(mp.nstr(mp.im(s), 55))
    JT = J.transpose()
    mr, mi = JT.solve(lr), JT.solve(li)
    tk = [-(mp.mpf(mr[k, 0].mid().str(40, radius=False)) + 1j * mp.mpf(mi[k, 0].mid().str(40, radius=False)))
          * mp.mpf(F[k, 0].mid().str(40, radius=False)) for k in range(N + 1)]
    c = mp.fsum(tk)
    cref = mp.mpc(*[mp.mpf(x) for x in am[str(N)]["c"]])
    absum = mp.fsum(abs(x) for x in tk)
    # Anteil der letzten Gleichungen
    teil = {K: mp.fsum(tk[N - K + 1:]) / c for K in (1, 2, 5, 10, 20, N // 2)}
    imax = max(range(N + 1), key=lambda k: abs(tk[k]))
    mabs = [abs(mp.mpf(mr[k, 0].mid().str(20, radius=False)) + 1j * mp.mpf(mi[k, 0].mid().str(20, radius=False))) for k in range(N + 1)]
    Fabs = [abs(mp.mpf(F[k, 0].mid().str(20, radius=False))) for k in range(N + 1)]
    print(f"N = {N}: c^lin = −⟨m, F⟩ = {mp.nstr(c, 8)}; c aus p_N (k7) = {mp.nstr(cref, 8)}; rel. Abstand {mp.nstr(abs(c - cref) / abs(cref), 3)}"
          f"   ({time.time() - t0:.0f} s)")
    print(f"   Auslöschung Σ|m_k F_k| / |c| = {mp.nstr(absum / abs(c), 3)}; grösster Beitrag bei k = {imax} ({mp.nstr(abs(tk[imax]) / abs(c), 3)}·|c|)")
    print("   Anteil der obersten K Gleichungen an c: " + ", ".join(f"K={K}: {mp.nstr(v, 4)}" for K, v in teil.items()))
    print("   |m_k|·|L|^{−k} bei k = 0, N/4, N/2, 3N/4, N: " + ", ".join(mp.nstr(mabs[k] * RL ** (-k), 3) for k in (0, N // 4, N // 2, 3 * N // 4, N))
          + ";  |F_k|·|L|^k: " + ", ".join(mp.nstr(Fabs[k] * RL ** k, 3) for k in (0, N // 4, N // 2, 3 * N // 4, N)), flush=True)
    # a_0-Jet-Basis (Strang B: J = H(S+Γ), F = −r = −H j): c = ⟨H^T m, j⟩ mit j = H^{-1}(−F)
    h = arb_poly([arb(0)] + a[1:])
    Hm = arb_mat(N + 1, N + 1)
    pw = arb_poly([1])
    for i in range(N + 1):
        for k in range(N + 1):
            Hm[k, i] = pw[k]
        pw = (pw * h).truncate(N + 1)
    jv = Hm.solve(F)
    mtr, mti = Hm.transpose() * mr, Hm.transpose() * mi
    uk = [-(mp.mpf(mtr[i, 0].mid().str(40, radius=False)) + 1j * mp.mpf(mti[i, 0].mid().str(40, radius=False)))
          * mp.mpf(jv[i, 0].mid().str(40, radius=False)) for i in range(N + 1)]
    cu = mp.fsum(uk)
    absu = mp.fsum(abs(x) for x in uk)
    imx = max(range(N + 1), key=lambda i: abs(uk[i]))
    jabs = [abs(mp.mpf(jv[i, 0].mid().str(20, radius=False))) for i in range(N + 1)]
    mt = [abs(mp.mpf(mtr[i, 0].mid().str(20, radius=False)) + 1j * mp.mpf(mti[i, 0].mid().str(20, radius=False))) for i in range(N + 1)]
    print(f"   a0-Jet-Basis: c = {mp.nstr(cu, 8)}; Auslöschung Σ|m~_i j_i|/|c| = {mp.nstr(absu / abs(c), 3)}; grösster Beitrag i = {imx}")
    for rho in (mp.mpf("0.5"), mp.mpf("0.876"), mp.mpf("1.349")):
        print(f"     ρ1 = {rho}: Σ|j_i|ρ1^i = {mp.nstr(mp.fsum(x * rho ** i for i, x in enumerate(jabs)), 3)}, max|m~_i|ρ1^−i = "
              f"{mp.nstr(max(x * rho ** (-i) for i, x in enumerate(mt)), 3)}", flush=True)
