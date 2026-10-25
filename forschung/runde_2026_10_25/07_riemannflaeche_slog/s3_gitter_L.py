"""s3: Log-Verzweigung an den Gitterpunkten L + 2 pi i k (Rand der periodischen Halbebene) im Motor-Blatt.
Test 1: k(w0+e1) - k(w0+e2) gegen (log e1 - log e2)/L  (Leitglied log(w-w0)/L).
Test 2: Funktionalgleichung k(w) = k(e^w) - 1 an den Kreispunkten (gleiches Blatt?).
Aufruf: python s3_gitter_L.py <dps>"""
import sys
import mpmath as mp
import motor

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = 40
L = mp.mpc("0.3181315052047641353126542515876645172035", "1.33723570143068940890116214319371061254")
pk = []
lab = []
for kk in [0, 1, 2]:
    w0 = L + 2j * mp.pi * kk
    for r in [mp.mpf("0.01"), mp.mpf("0.001")]:
        for phi in [mp.pi, mp.pi / 2, -mp.pi / 2, mp.mpf(0)]:
            w = w0 + r * mp.expjpi(phi / mp.pi)
            pk.append(w)
            lab.append((kk, r, phi))
aus = [f"slog({motor.gpzahl(w, 30)})" for w in pk] + [f"slog({motor.gpzahl(mp.exp(w), 30)})" for w in pk]
res = motor.werte("exp(1)", aus, dps)
n = len(pk)
kv = res[:n]
ke = res[n:]
for i, (kk, r, phi) in enumerate(lab):
    print(f"k={kk} r={mp.nstr(r,2)} phi={mp.nstr(phi,4):>6}: slog={mp.nstr(kv[i],15):>40}  FG-Rest k(w)-k(e^w)+1={mp.nstr(abs(kv[i]-ke[i]+1),3)}")
print("Log-Test: Differenz r=0.001 minus r=0.01 bei gleichem Winkel; Theorie log(0.1)/L =", mp.nstr(mp.log(mp.mpf('0.1')) / L, 10))
for i in range(n):
    kk, r, phi = lab[i]
    if r == mp.mpf("0.01"):
        j = i + 4
        print(f"  k={kk} phi={mp.nstr(phi,4):>6}: {mp.nstr(kv[j]-kv[i],10)}   Abw. {mp.nstr(abs(kv[j]-kv[i]-mp.log(mp.mpf('0.1'))/L),3)}")
motor.schliessen()
