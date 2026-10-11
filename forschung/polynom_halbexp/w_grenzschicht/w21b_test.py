import sys, numpy as np
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", "N = 300"))
for z in [3.0, 5.0, 10.0, -3.0, 2.0+2.0j]:
    v = T(z)
    print(f"ζ={z}: (1/N)ln|T| = {v:.4f}, −log|ζ|·(N+2)/N = {-np.log(abs(z))*(N+2)/N:.4f}")
