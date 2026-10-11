"""Reeller Strahl eines Log-Punkts Z von F_q: von F_q = −6 (nahe Z) per Newton-Fortsetzung F_q(ζ) = t, t ↑ L_T. Schliesst er bei Z + 1?
 (Kneser-Kontrolle: Z = ζ_c − 2 muss bei ζ_c − 1 ankommen.) F_q über die Hauptkette des Forks (nur gelesen).
 Aufruf: strahl.py "Re1,Im1;Re2,Im2;..."
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
from fatou_backend import FatouGP, find_default_gp_exe

LT = -1.743909176132737
t0 = time.time()
punkte = [complex(*map(float, s.split(","))) for s in sys.argv[1].split(";")]
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)


def newton(z, w):
    v = gp.eval_batch("2", [f"if(1, kr = iferr(if(1, kv=alarm(4, kNewt({z.real:.20f}+I*({z.imag:.20f}), {w})); if(type(kv)==\"t_ERROR\", [10^9, 10^9], kv)),E,[10^9, 10^9]); kr[1])", "kr[2]"])
    return complex(v[0]), abs(complex(v[1]))


try:
    gp.eval_batch("2", ["hexinit(sexp(1),0,1)",
                        "if(1, kpsi = Pol(serreverse(ksig + O(x^(poldegree(ksig)+1)))); kfrs = abs(subst(ksig, x, krs))/2; 0)",
                        "if(1, kFq = ((z) -> if(1, my(s = exp(z*klog), n = 0, y); while(abs(s) > kfrs, s = s/klam; n++); y = kp + subst(kpsi, x, s); for(i = 1, n, y = pent(y)); y)); 0)",
                        "if(1, kNewt = ((z, w) -> if(1, my(f, d); for(i = 1, 25, f = kFq(z) - w; d = (kFq(z + 1e-15) - kFq(z - 1e-15))/2e-15; z = z - f/d; if(abs(f) < 1e-25, break)); [z, abs(kFq(z) - w)])); 0)"])
    for Z in punkte:
        # Startpunkt auf dem Strahl: Vorzeichenwechsel von Im F_q auf dem Kreis |ζ − Z| = 1e−3 (stetiger Teil, nicht der Schnitt der Hauptkette)
        phis = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        kreis = Z + 1e-3 * np.exp(1j * phis)
        F = np.array([complex(v) for v in gp.eval_batch("2", [f"kFq({u.real:.18f}+I*({u.imag:.18f}))" for u in kreis])])
        wechsel = [i for i in range(64) if F[i].imag * F[(i + 1) % 64].imag < 0 and abs(F[i] - F[(i + 1) % 64]) < 1.0]
        if len(wechsel) != 1:
            print(f"Z = {Z:.6f}: {len(wechsel)} Vorzeichenwechsel auf dem Kreis — kein eindeutiger Strahl [{time.time() - t0:.0f}s]", flush=True)
            continue
        i = wechsel[0]
        t_start = float(F[i].real)
        z, r = newton(kreis[i], t_start)
        if r > 1e-12:
            print(f"Z = {Z:.6f}: Newton am Strahlanfang gescheitert [{time.time() - t0:.0f}s]", flush=True)
            continue
        weg, smax = [z], 0.0
        ts = np.concatenate([np.linspace(t_start, -2.2, 90), np.linspace(-2.2, LT - 1e-4, 60)])
        for t in ts[1:]:
            zn, r = newton(z, t)
            if r > 1e-12:
                print(f"   Newton bricht ab bei t = {t:.4f} (ζ = {z:.6f})", flush=True)
                break
            smax = max(smax, abs(zn - z))
            z = zn
            weg.append(z)
        print(f"Z = {Z:.6f}: Strahl von {weg[0]:.6f} bis {z:.6f} (t = {t:.4f}), grösster Schritt {smax:.3f}; "
              f"Abstand zu Z + 1: {abs(z - (Z + 1)):.4f} [{time.time() - t0:.0f}s]", flush=True)
        np.save(f"strahl_{Z.real:+.3f}_{Z.imag:+.3f}.npy", np.array(weg))
finally:
    gp.close()
