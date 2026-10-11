"""Gebietsunabhängigkeit an einer Familie von Fundamentalgebieten G_c: γ_c = Kreisbogen von q̄ über den reellen Punkt c nach q, G_c zwischen γ_c und P(γ_c).
 Für jedes c: (a) γ_c ∩ P(γ_c) = {q, q̄}? (b) Γ_c = α_K(γ_c) (Kneser-Abel-Funktion; unten ihex mit Zeitgrenze, oben θ-Darstellung des Forks);
 enthält der Streifen zwischen Γ_c und Γ_c + 1 einen Verzweigungspunkt z_{m,k} = z₀ + m + (k − 1)Π von H_K (z₀ = −0.7408 + 2.7458i)?
 Wenn (a) und nicht (b): H_K bildet den Streifen biholomorph auf G_c ab, α_K ist auf G_c holomorph und injektiv → Eindeutigkeit auch auf G_c.
 Fork nur gelesen. Aufruf: gebiete_gc.py c1,c2,...
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

cs = [mp.mpf(c) for c in sys.argv[1].split(",")]
dps = 30
mp.mp.dps = dps + 10
t0 = time.time()
PI_ = complex(3.00910777634078, 2.9413168766409)
Z0 = complex(-0.7407566463, 2.745812873)
VERZW = [Z0 + m + (k - 1) * PI_ for k in (1, 2, 3) for m in range(0, 12)]
VERZW += [v.conjugate() for v in VERZW]
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)


def bogen(c, q, n):
    """Punkte auf dem Kreisbogen von c (t = 0) nach q (t = 1), symmetrisch zur reellen Achse (Mittelpunkt auf ℝ)."""
    if abs(c - q.real) < 1e-12:
        return [complex(c, q.imag * t) for t in np.linspace(0, 1, n)]
    m = (c * c - abs(q) ** 2) / (2 * (c - q.real))
    r = abs(c - m)
    a0, a1 = np.angle(c - m), np.angle(q - m)
    if a1 - a0 > np.pi:
        a1 -= 2 * np.pi
    if a0 - a1 > np.pi:
        a1 += 2 * np.pi
    return [m + r * np.exp(1j * (a0 + (a1 - a0) * t)) for t in np.linspace(0, 1, n)]


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


def im_streifen(G, z):
    """Liegt z im Streifen zwischen der Kurve G (von der Achse nach oben) und G + 1? Waagrechter Schnitt auf Höhe Im z."""
    for i in range(len(G) - 1):
        a, b = G[i], G[i + 1]
        if (a.imag - z.imag) * (b.imag - z.imag) <= 0 and a.imag != b.imag:
            x = a.real + (z.imag - a.imag) * (b.real - a.real) / (b.imag - a.imag)
            if x < z.real < x + 1:
                return True
    return False


try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", "real(kp)", "imag(kp)"])
    q = complex(k[1].real, k[2].real)
    print(f"q = {q:.6f}, Stellen {mp.nstr(k[0], 4)} [{time.time() - t0:.0f}s]", flush=True)
    for c in cs:
        pts = bogen(float(c), q, 160)[:-1]
        bild = [complex(v) for v in gp.eval_batch("2", [f"iferr(pent({w.real:.25f}+I*{w.imag:.25f}),E,10^9)" for w in pts])]
        G = np.array(pts); PG = np.array(bild)
        Gv = np.concatenate([np.conj(G[::-1]), G[1:]]); PGv = np.concatenate([np.conj(PG[::-1]), PG[1:]])
        s_a = schnitte(Gv, PGv)
        # α_K: oben θ-Darstellung, unten ihex (Zeitgrenze 2 s), Wechsel wo beide da sind und übereinstimmen
        liste = "[" + ",".join(f"{w.real:.25f}+I*{w.imag:.25f}" for w in pts) + "]"
        ob = gp.eval_batch("2", [f"if(1, kgp={liste}; kgW=knwalk(kgp); kgA=vector(#kgp,j,kgW[j]+sum(n=1,#kncc,kncc[n]*knen(kgW[j],n-1))); 0)"]
                           + [f"kgA[{j + 1}]" for j in range(len(pts))] + [f"kns*imag(kgW[{j + 1}])" for j in range(len(pts))])
        A_th = [complex(v) for v in ob[1:len(pts) + 1]]
        h = [v.real for v in ob[len(pts) + 1:]]
        A_ih = [complex(v) for v in gp.eval_batch("2", [f"iferr(if(1, kv=alarm(2, ihex({w.real:.25f}+I*{w.imag:.25f})); if(type(kv)==\"t_ERROR\", 10^9, kv)),E,10^9)" for w in pts])]
        GA, wechsel, abw = [], None, None
        for j in range(len(pts)):
            if abs(A_ih[j]) < 1e8 and (wechsel is None):
                GA.append(A_ih[j])
                if h[j] > 1.9 and abs(A_ih[j] - A_th[j]) < 1e-10:
                    wechsel, abw = j, abs(A_ih[j] - A_th[j])
            else:
                if wechsel is None:
                    wechsel, abw = j, (abs(A_ih[j - 1] - A_th[j - 1]) if j else None)
                GA.append(A_th[j])
        GA = np.array(GA)
        drin = [v for v in VERZW if im_streifen(GA, v)]
        ok = (s_a == 0) and not drin and abs(GA[0].imag) < 1e-8
        print(f"c = {mp.nstr(c, 5):7s}: (a) Schnitte {s_a}; Wechsel ihex→θ bei Punkt {wechsel} (Abweichung {abw if abw is None else f'{abw:.1e}'}); "
              f"α_K(c) = {GA[0].real:.5f}; max Re Γ_c = {GA.real.max():.4f}; Verzweigungspunkte im Streifen: "
              f"{[f'{v:.3f}' for v in drin] if drin else 'keine'} → {'α_K holomorph und injektiv auf G_c' if ok else 'nicht gezeigt'} [{time.time() - t0:.0f}s]", flush=True)
finally:
    gp.close()
