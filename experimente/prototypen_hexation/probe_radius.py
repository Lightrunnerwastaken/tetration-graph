"""Konvergenzradius der Taylorreihe A der Abel-Funktion α (Umkehrung der Hexation) um c = Re q im Fork; Ort der nächsten Singularität.
 a_k reell; bei einem konjugierten Paar w* = c + R e^{±iφ} gilt a_{k+1} ≈ (2 cos φ / R) a_k − a_{k−1}/R² — Fit über das Ende der Reihe."""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis = sys.argv[1] if len(sys.argv) > 1 else "2"
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = dps + 10
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
try:
    kopf = gp.eval_batch(basis, ["hexinit(sexp(1),0,1)", "knc", "knrho", "poldegree(knA)", "hexerr", "knr"])
    c, rho, grad, err, knr = kopf[1].real, kopf[2].real, int(kopf[3].real), kopf[4].real, kopf[5].real
    a = [w.real for w in gp.eval_batch(basis, ["hexinit(sexp(1),0,1)"] + [f"polcoeff(knA,{k})" for k in range(grad + 1)])[1:]]
finally:
    gp.close()
print(f"Basis {basis}: c = {mp.nstr(c, 12)}, ρ = {mp.nstr(rho, 8)}, Grad {grad}, Stellen {mp.nstr(err, 5)}, Im p = {mp.nstr(knr, 8)}")
b = [a[k] / rho ** k for k in range(grad + 1)]            # Koeffizienten in (w − c)
for k in range(0, grad + 1, max(1, grad // 12)):
    print(f"  k = {k:4d}: |b_k| = {mp.nstr(abs(b[k]), 6)},  |b_k|^(-1/k) = {mp.nstr(abs(b[k]) ** (-mp.mpf(1) / k), 6) if k and b[k] else '—'}")
# Zwei-Glieder-Rekursion b_{k+1} = s b_k + t b_{k−1} per kleinste Quadrate über ein Fenster
for lo, hi in ((grad // 4, grad // 2), (grad // 2, 3 * grad // 4), (grad // 3, 2 * grad // 3)):
    A = mp.matrix([[b[k], b[k - 1]] for k in range(lo, hi)])
    y = mp.matrix([b[k + 1] for k in range(lo, hi)])
    s, t = mp.lu_solve(A.T * A, A.T * y)
    # z² = s z + t → Wurzeln 1/w-Richtung: b_k ~ r^k mit r = Wurzel; Singularität bei 1/r
    disk = s * s + 4 * t
    r1 = (s + mp.sqrt(disk)) / 2
    w = c + 1 / r1
    print(f"  Fenster k = {lo}..{hi}: Singularität bei w* ≈ {mp.nstr(w, 8)}, |w* − c| = {mp.nstr(abs(1 / r1), 6)} = {mp.nstr(abs(1 / r1) / knr, 4)} Im p")
