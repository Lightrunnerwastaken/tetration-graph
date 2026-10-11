"""Runde 2, Schritt 7: Brouwer-Reduktion bei exp prüfen. Sekanten-Operator 𝒥_p u = u∘p + M_p·u,
M_p = (f∘p − f∘f)/(p − f) = Σ_k f_k h_k, h_k = Σ_{j<k} p^j f^{k−1−j} (h_{k+1} = p h_k + f^k)  [keine Division durch p − f ≈ 0].
T(p) := Π^{(p)} f_K = (T_N𝒥_p|_{P_N})⁻¹ T_N[f_K∘p + M_p f_K].  Für p = p_N gilt T(p_N) = p_N (Fixpunkt).
Frage: bleibt ‖T(p) − f_K‖ ≈ ‖p_N − f_K‖ (winzig), wenn p in einer festen Kugel um p_N variiert?  Dann trägt der Brouwer-Schluss.
Daten: gemeinsam/kneser_halbexp_koeff_1000.json (Kneser, ~47 Stellen), gemeinsam/polynome_bis_194.json (p_N)."""
import json, sys, time, random
from pathlib import Path
import mpmath as mp
from flint import arb, acb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).parent
G = HIER.parent / "gemeinsam"
mp.mp.dps = 50
ctx.prec = 700
Ns = [int(x) for x in sys.argv[1].split(",")]
etas = [mp.mpf(x) for x in sys.argv[2].split(",")]
KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 700
fk = [arb(x) for x in json.loads((G / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
PN = json.loads((G / "polynome_bis_194.json").read_text())
PN.update(json.loads((G / "polynome_195_bis_300.json").read_text()))
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
rL = abs(L)
Kz = 256


def fK_wert(z):
    s = mp.mpc(0)
    for k in range(len(fk) - 1, -1, -1):
        s = s * z + mp.mpf(fk[k].mid().str(50, radius=False))
    return s


def T_von_p(pk, N):
    p = arb_poly(pk)
    f = arb_poly(fk[:N + 1])
    P, h, F, M, fop = arb_poly([1]), arb_poly([1]), arb_poly([1]), arb_poly([0]), arb_poly([fk[0]])
    Pj = [arb_poly([1])]
    for k in range(1, KMAX):
        # h = h_k, F = f^{k−1}, P = p^{k−1}
        M += fk[k] * h
        P = (P * p).truncate(N + 1)                 # p^k
        fop += fk[k] * P
        if k <= N:
            Pj.append(P)
        F = (F * f).truncate(N + 1)                 # f^k
        h = ((p * h).truncate(N + 1) + F).truncate(N + 1)   # h_{k+1}
    M = M.truncate(N + 1)
    A = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = Pj[j] + (M * arb_poly([0] * j + [1])).truncate(N + 1)
        for k in range(N + 1):
            A[k, j] = col[k]
    rhs = (fop + (M * f).truncate(N + 1)).truncate(N + 1)
    b = arb_mat([[rhs[k]] for k in range(N + 1)])
    x = A.solve(b)
    return [mp.mpf(x[k, 0].mid().str(48, radius=False)) for k in range(N + 1)]


random.seed(7)
for N in Ns:
    t0 = time.time()
    pN = [mp.mpf(x) for x in PN[str(N)]]
    zs = [rL * mp.mpf("0.999") * mp.expjpi(2 * mp.mpf(j) / Kz) for j in range(Kz)]
    fz = [fK_wert(z) for z in zs]
    err_pN = max(abs(mp.polyval(pN[::-1], z) - v) for z, v in zip(zs, fz))
    print(f"N={N}: max_(|z|=0.999|L|)|p_N − f_K| = {mp.nstr(err_pN, 3)}", flush=True)
    for eta in etas:
        q = [random.gauss(0, 1) * float(rL) ** (-k) / (k + 1) ** 2 for k in range(N + 1)]
        qn = sum(abs(q[k]) * float(rL) ** k for k in range(N + 1))
        pk = [arb(mp.nstr(pN[k] + eta * q[k] / qn, 55)) for k in range(N + 1)]
        Tp = T_von_p(pk, N)
        err = max(abs(mp.polyval(Tp[::-1], z) - v) for z, v in zip(zs, fz))
        dfix = max(abs(mp.polyval(Tp[::-1], z) - mp.polyval(pN[::-1], z)) for z in zs)
        print(f"   ‖p − p_N‖ ≤ η = {mp.nstr(eta, 2)} (Gewicht |L|^k):  max|T(p) − f_K| = {mp.nstr(err, 3)},  max|T(p) − p_N| = {mp.nstr(dfix, 3)}"
              f"  ({time.time() - t0:.0f}s)", flush=True)
