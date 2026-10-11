"""Fortsetzung der Kneser-Abel-Funktion α_K = H_K⁻¹ längs eines Weges: reeller Punkt x0 → hoch zum Kreis um p → Umlauf → zurück nach x0.
 Newton auf hex(z) = w(t) (Fork), Ableitung über zentrale Differenzen; je Schritt Kontrolle |hex(z) − w| und Stetigkeit.
 Ergebnis α_K^neu(x0): reell → keine neue Lösung; komplex → H_K ist dort reell (neue reelle Kurve oder Kopie).
 Fork nur gelesen. Aufruf: fortsetzung_p.py x0 r [umläufe]
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

x0, r = float(sys.argv[1]), float(sys.argv[2])
umlaeufe = int(sys.argv[3]) if len(sys.argv) > 3 else 1
mp.mp.dps = 40
t0 = time.time()
P_FIX = complex(1.7937303432391165676, 1.5736302435118265996)
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
H = 1e-12


def hexv(zs):
    ex = []
    for z in zs:
        for d in (0, H, -H):
            ex.append(f"iferr(if(1, kv=alarm(3, hex({z.real + d:.25f}+I*{z.imag:.25f})); if(type(kv)==\"t_ERROR\", 10^9, kv)),E,10^9)")
    v = [complex(x) for x in gp.eval_batch("2", ex)]
    return [(v[3 * i], (v[3 * i + 1] - v[3 * i + 2]) / (2 * H)) for i in range(len(zs))]


def folge(z, weg, name):
    """z = α_K(weg[0]); Newton-Fortsetzung längs weg. Rückgabe Endwert und Protokoll."""
    zs, maxres = [z], 0.0
    for j, w in enumerate(weg[1:], 1):
        zn = z
        for _ in range(6):
            (f, df), = hexv([zn])
            if abs(f) > 1e8 or df == 0:
                raise RuntimeError(f"{name}: hex versagt bei z = {zn:.6f} (Schritt {j})")
            dz = (f - w) / df
            zn = zn - dz
            if abs(dz) < 1e-13:
                break
        (f, _), = hexv([zn])
        maxres = max(maxres, abs(f - w))
        if abs(zn - z) > 0.3:
            print(f"   Warnung {name}: Sprung {abs(zn - z):.3f} bei Schritt {j} (w = {w:.4f}, z {z:.4f} → {zn:.4f})", flush=True)
        z = zn
        zs.append(z)
    print(f"  {name}: {len(weg)} Punkte, grösster Rest |hex(z) − w| = {maxres:.1e}, Ende z = {z:.10f} [{time.time() - t0:.0f}s]", flush=True)
    return z, np.array(zs)


try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", f"ihex({x0})"])
    z = complex(k[1])
    print(f"x0 = {x0}: α_K(x0) = {z:.12f} (reell) [{time.time() - t0:.0f}s]", flush=True)
    w_start = P_FIX - r                                        # linker Kreispunkt (φ = 180°)
    hoch = [x0 + (w_start - x0) * t for t in np.linspace(0, 1, 80)]
    z, _ = folge(z, hoch, "hoch")
    z_kreis0 = z
    phis = np.pi + np.linspace(0, 2 * np.pi * umlaeufe, 144 * umlaeufe + 1)
    kreis = [P_FIX + r * np.exp(1j * f) for f in phis]
    z, zk = folge(z, kreis, f"{umlaeufe}× um p (gegen den Uhrzeigersinn)")
    print(f"  Monodromie am Kreispunkt: Δ = {z - z_kreis0:.10f}", flush=True)
    runter = hoch[::-1]
    z, _ = folge(z, runter, "zurück")
    print(f"ERGEBNIS: α_K^neu({x0}) = {z:.12f}  (vorher {complex(k[1]).real:.12f})", flush=True)
    chk = gp.eval_batch("2", [f"hex({z.real:.25f}+I*{z.imag:.25f})"])[0]
    print(f"  Probe hex(α_K^neu) = {mp.nstr(chk, 15)}", flush=True)
finally:
    gp.close()
