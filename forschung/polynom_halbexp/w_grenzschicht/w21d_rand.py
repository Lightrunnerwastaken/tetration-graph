import sys, numpy as np
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", "N = 300"))
a0 = coef[0]; R0 = abs(L - a0); psiL = np.angle(L - a0); argL = np.angle(L)
vals = []
for p in np.linspace(psiL, -psiL, 61):
    vals.append(("A", np.degrees(p), T(a0 + R0*np.exp(1j*p))))
for ph in np.linspace(argL, 2*np.pi - argL, 61):
    vals.append(("D", np.degrees(ph), T(rho*np.exp(1j*ph))))
mx = max(vals, key=lambda v: v[2])
print("max (1/N)ln|T| auf ∂W:", mx)
for v in vals[::6]: print(v)
