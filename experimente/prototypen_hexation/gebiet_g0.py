"""Festes Fundamentalgebiet G₀ für die Eindeutigkeit (Trappmann–Kouznetsov-Mechanismus ohne Gebietslücke).
 γ₀ = senkrechte Strecke von q̄ nach q durch Re q; G₀ = Gebiet zwischen γ₀ und P(γ₀).
 Prüft (a) γ₀ ∩ P(γ₀) = {q, q̄} (Jordan-Rand), (b) Γ₀ = α_K(γ₀) (Kneser-Abel-Funktion) hat Re < −1.741, so dass der Streifen zwischen
 Γ₀ und Γ₀ + 1 im singularitätenfreien Gebiet Re z < −0.741 liegt und H_K ihn biholomorph auf G₀ abbildet.
 α_K unten über ihex (Pentation bis in den Bogen), oben über die θ-Darstellung des Forks: α = W + Σ cₙ e^{2πi s (n−1) W}, W = α_p (knwalk).
 Fork nur gelesen. Aufruf: gebiet_g0.py [dps]
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fatou_backend import FatouGP, find_default_gp_exe

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
N_U, N_O = 31, 140                       # Punkte unten (|t| ≤ 0.3 Im q, ihex) und oben (0.25 … 0.995 Im q, θ-Darstellung)
try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", "real(kp)", "imag(kp)", "hexerr"])
    req, imq = k[1].real, k[2].real
    print(f"Stellen {mp.nstr(k[3], 5)}, Re q = {mp.nstr(req, 10)}, Im q = {mp.nstr(imq, 8)}", flush=True)
    t_u = [imq * mp.mpf("0.3") * mp.mpf(j) / (N_U - 1) for j in range(N_U)]
    unten = gp.eval_batch("2", [f"ihex({mp.nstr(req, dps)}+I*{mp.nstr(t, dps)})" for t in t_u])
    print(f"  unten (ihex): {N_U} Punkte [{time.time() - t0:.0f}s]", flush=True)
    t_o = [imq * (mp.mpf("0.25") + (mp.mpf("0.995") - mp.mpf("0.25")) * mp.mpf(j) / (N_O - 1)) for j in range(N_O)]
    liste = "[" + ",".join(f"{mp.nstr(req, dps)}+I*{mp.nstr(t, dps)}" for t in t_o) + "]"
    oben = gp.eval_batch("2", [f"if(1, kgp={liste}; kgW=knwalk(kgp); "
                                "kgA=vector(#kgp,j,kgW[j]+sum(n=1,#kncc,kncc[n]*knen(kgW[j],n-1))); 0)"]
                         + [f"kgA[{j + 1}]" for j in range(N_O)] + [f"kns*imag(kgW[{j + 1}])" for j in range(N_O)])
    A_o = [v for v in oben[1:N_O + 1]]
    h_o = [v.real for v in oben[N_O + 1:]]
    print(f"  oben (θ-Darstellung): s·Im W von {mp.nstr(min(h_o), 4)} bis {mp.nstr(max(h_o), 4)} [{time.time() - t0:.0f}s]", flush=True)
    for f in ("0.25", "0.27", "0.3"):
        t = imq * mp.mpf(f)
        j = min(range(N_O), key=lambda i: abs(t_o[i] - t))
        a = gp.eval_batch("2", [f"ihex({mp.nstr(req, dps)}+I*{mp.nstr(t_o[j], dps)})"])[0]
        print(f"  Abgleich t = {mp.nstr(t_o[j] / imq, 4)} Im q: ihex {mp.nstr(a, 12)}, θ-Darstellung {mp.nstr(A_o[j], 12)}, Unterschied {mp.nstr(abs(a - A_o[j]), 3)}", flush=True)
    t_p = [imq * mp.mpf(j) / 400 for j in range(0, 399)]
    bild = gp.eval_batch("2", [f"pent({mp.nstr(req, dps)}+I*{mp.nstr(t, dps)})" for t in t_p])
finally:
    gp.close()
G0 = np.array([complex(req, float(t)) for t in t_p])
PG = np.array([complex(b) for b in bild])
GA = np.array([complex(a) for a in unten] + [complex(a) for a, t in zip(A_o, t_o) if t > t_u[-1]])


def schnitte(c1, c2):
    n = 0
    for i in range(len(c1) - 1):
        for j in range(len(c2) - 1):
            p1, p2, p3, p4 = c1[i], c1[i + 1], c2[j], c2[j + 1]
            d = (p2 - p1).real * (p4 - p3).imag - (p2 - p1).imag * (p4 - p3).real
            if d == 0:
                continue
            t = ((p3 - p1).real * (p4 - p3).imag - (p3 - p1).imag * (p4 - p3).real) / d
            s = ((p3 - p1).real * (p2 - p1).imag - (p3 - p1).imag * (p2 - p1).real) / d
            if 1e-12 < t < 1 - 1e-12 and 1e-12 < s < 1 - 1e-12:
                n += 1
    return n


G0v = np.concatenate([np.conj(G0[::-1]), G0[1:]])
PGv = np.concatenate([np.conj(PG[::-1]), PG[1:]])
print(f"(a) Schnitte γ₀ mit P(γ₀) (obere und untere Hälfte, ohne Enden): {schnitte(G0v, PGv)}; P(γ₀) Selbstschnitte: "
      f"{sum(schnitte(PGv[i:i+2], PGv[j:j+2]) for i in range(len(PGv) - 1) for j in range(i + 2, len(PGv) - 1))}")
print(f"    P(γ₀) oben: von {mp.nstr(PG[0], 8)} nach {mp.nstr(PG[-1], 8)} (q = {mp.nstr(mp.mpc(req, imq), 8)})")
print(f"(b) max Re Γ₀ = {GA.real.max():.6f} (muss < −1.741); Re Γ₀ oben {GA[-1].real:.4f}, Im Γ₀ oben {GA[-1].imag:.4f}")
fig, ax = plt.subplots(1, 2, figsize=(13, 6))
ax[0].plot(G0v.real, G0v.imag, color="#3d6fb6", label="γ₀ (senkrecht durch Re q)")
ax[0].plot(PGv.real, PGv.imag, color="#b5523b", label="P(γ₀)")
for z, m in ((complex(req, imq), "ko"), (complex(req, -imq), "ko"), (complex(-1.9786, 1.91185), "k^")):
    ax[0].plot(z.real, z.imag, m)
ax[0].set_aspect("equal"); ax[0].legend(fontsize=8); ax[0].set_title("G₀ zwischen γ₀ und P(γ₀) (w-Ebene)")
ax[1].plot(GA.real, GA.imag, color="#3d6fb6", label="Γ₀ = α_K(γ₀), obere Hälfte")
ax[1].plot(GA.real + 1, GA.imag, color="#b5523b", label="Γ₀ + 1")
ax[1].axvline(-0.741, ls=":", color="k", label="erste Verzweigungen ab Re z = −0.741")
ax[1].axvline(-3.767, ls="--", color="#8a93a3", label="s_b")
ax[1].set_aspect("equal"); ax[1].legend(fontsize=8); ax[1].set_title("Streifen Σ₀ = α_K(G₀) (z-Ebene)")
plt.tight_layout()
plt.savefig(r"D:/Users/janis/Documents/Matura/tetration-graph/experimente/prototypen_hexation/gebiet_g0.png", dpi=90)
print(f"fertig [{time.time() - t0:.0f}s]")
