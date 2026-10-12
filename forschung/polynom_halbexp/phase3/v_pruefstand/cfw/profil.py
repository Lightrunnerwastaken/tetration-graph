import numpy as np, math, sys
sys.stdout.reconfigure(encoding="utf-8")
datei = sys.argv[1]
W = np.load(datei)
t = W["t"]; L = complex(0.31813150520476413, 1.3372357014306895); rho = abs(L); thL = math.atan2(L.imag, L.real)
ang = np.mod(np.angle(t), 2*np.pi)
gin = (ang > thL) & (ang < 2*np.pi - thL)
dL = np.degrees(np.minimum(np.abs(np.angle(t/L)), np.abs(np.angle(t/np.conj(L)))))
dist = np.minimum(np.abs(t - L), np.abs(t - np.conj(L)))
s = np.log(np.abs(t)/rho)
r = {}
for n1, n2 in [(1000,1500),(1500,2000)]:
    r[(n1,n2)] = (W[f"lw_{n2}"] - W[f"lw_{n1}"])/(n2-n1)
print("Abstand zu L, L̄ (Grad) | |t−L| | Linsenhöhe s | Rate 1000→1500 | 1500→2000 | ln w(2000)")
for dmin, dmax in [(0.5,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,10),(10,15),(15,30),(30,60),(60,104)]:
    j = np.where(gin & (dL >= dmin) & (dL < dmax))[0]
    jm = j[np.argmax(r[(1500,2000)][j])]
    print(f"{dmin:>4}–{dmax:<4}° | {dist[jm]:.3f} | {s[jm]:.5f} | {r[(1000,1500)][jm]:+.5f} | {r[(1500,2000)][jm]:+.5f} | {W['lw_2000'][jm]:+.3f}  (Max im Band bei {math.degrees(ang[jm]):.1f}°)")
for key in [(1000,1500),(1500,2000)]:
    pos = gin & (r[key] > 0)
    if pos.any():
        print(key, f"Rate > 0 auf der Linse bis {dL[pos].max():.2f}° von L bzw. L̄ (|t − L| ≤ {dist[pos].max():.3f}, Linsenhöhe dort ≤ {s[pos].max():.5f})")
