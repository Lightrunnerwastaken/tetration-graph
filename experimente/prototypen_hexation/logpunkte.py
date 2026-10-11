"""Log-Punkte von F_q im Band zwischen C und C + Π: Z = ζ* + 1 mit F_q(ζ*) = w_b + 2πi·j/ln λ_T (w_b: Log-Punkt der Pentation, P → −∞).
 Jede reelle Lösung der θ-Klasse beginnt an so einem Z (einziger reeller Strahl); die Kneser-artige gehört zu Z = ζ_c − 2.
 F_q über die Hauptkette des Forks (nur gelesen; in karte_reell-Bereich gleich der verfolgten Fortsetzung, 1e−14).
 Ausgabe: Liste der Z und ob F_q dort stetig ist (Nachbarn) — Aufruf: logpunkte.py
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
from fatou_backend import FatouGP, find_default_gp_exe

t0 = time.time()
ZETA_C = complex(0.4031293, -1.7097883)
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
try:
    gp.eval_batch("2", ["hexinit(sexp(1),0,1)",
                        "if(1, kpsi = Pol(serreverse(ksig + O(x^(poldegree(ksig)+1)))); kfrs = abs(subst(ksig, x, krs))/2; 0)",
                        "if(1, kFq = ((z) -> if(1, my(s = exp(z*klog), n = 0, y); while(abs(s) > kfrs, s = s/klam; n++); y = kp + subst(kpsi, x, s); for(i = 1, n, y = pent(y)); y)); 0)",
                        # w_b: pent(w - 1) = -2 auf der Halbperiode
                        "if(1, kwb = -1.9786 + 1.91185*I; for(i = 1, 40, kwb = kwb - (pent(kwb - 1) + 2)/(pent(kwb - 1 + 1e-20) - pent(kwb - 1 - 1e-20))*2e-20); 0)",
                        "if(1, kNewt = ((z, w) -> if(1, my(f, d); for(i = 1, 15, f = kFq(z) - w; d = (kFq(z + 1e-15) - kFq(z - 1e-15))/2e-15; z = z - f/d; if(abs(f) < 1e-25, break)); [z, abs(kFq(z) - w)])); 0)"])
    wb = complex(gp.eval_batch("2", ["kwb"])[0])
    per = complex(gp.eval_batch("2", ["2*Pi*I/log(lambdat)"])[0])
    print(f"w_b = {wb:.15f}, Periode {per:.6f}; Probe F_q(ζ_c − 3) = {complex(gp.eval_batch('2', [f'kFq({ZETA_C.real - 3}+I*({ZETA_C.imag}))'])[0]):.6f} [{time.time() - t0:.0f}s]", flush=True)
    funde = []
    starts = [complex(x, y) for x in np.arange(-3.0, 6.01, 0.5) for y in np.arange(-1.6, 1.21, 0.4)]
    for j in (0, 1, -1, 2):
        ziel = wb + j * per
        ex = []
        for s in starts:
            ex += [f"if(1, kr = iferr(if(1, kv=alarm(4, kNewt({s.real}+I*({s.imag}), {ziel.real:.25f}+I*({ziel.imag:.25f}))); if(type(kv)==\"t_ERROR\", [10^9, 10^9], kv)),E,[10^9, 10^9]); kr[1])", "kr[2]"]
        v = [complex(x) for x in gp.eval_batch("2", ex)]
        for z, r in zip(v[0::2], v[1::2]):
            if abs(r) < 1e-12 and -3.2 <= z.real <= 6.5 and -1.75 <= z.imag <= 1.25 and all(abs(z - f[0]) > 1e-6 for f in funde):
                funde.append((z, j))
        print(f"  j = {j:+d}: bisher {len(funde)} Lösungen [{time.time() - t0:.0f}s]", flush=True)
finally:
    gp.close()
funde.sort(key=lambda t: (t[0].imag, t[0].real))
print("Log-Punkte Z = ζ* + 1 (Band zwischen C und C + Π, Hauptkette):")
for z, j in funde:
    Z = z + 1
    k_kneser = min(range(-2, 4), key=lambda k: abs(Z - (ZETA_C - 2 + k * (3.00910778 + 2.94131688j))))
    print(f"   Z = {Z:.8f}  (j = {j:+d}){'   = ζ_c − 2 + kΠ' if abs(Z - (ZETA_C - 2 + k_kneser * (3.00910778 + 2.94131688j))) < 1e-4 else ''}")
print(f"fertig [{time.time() - t0:.0f}s]")
