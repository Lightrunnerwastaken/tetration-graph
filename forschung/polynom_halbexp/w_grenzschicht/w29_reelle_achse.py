import sys, numpy as np
N = int(sys.argv[1])
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", f"N = {N}"))
a0 = coef[0]; R0 = abs(L - a0)
print(f"A-Scheitel bei x = {a0+R0:.4f} = {(a0+R0)/rho:.4f} ρ")
xs = np.array([1.02, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.33, 1.345, 1.36, 1.4, 1.5, 1.7, 2.0, 2.5]) * rho
for x in xs:
    print(f"  x/ρ = {x/rho:.3f}: (1/N)ln|T| = {T(complex(x, 1e-9)):+.4f}", flush=True)
