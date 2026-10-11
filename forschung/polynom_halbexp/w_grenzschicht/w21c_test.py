import sys, numpy as np
exec(open("w21_eckanteil_gitter.py", encoding="utf-8").read().split("radii = ")[0].replace("N = int(sys.argv[1])", "N = 150").replace("ctx.cap = N + 10", "ctx.cap = N + 10"))
# Residuum J^T Z − D u
res = J.transpose() * Z
r = max(abs(float((res[k,0] - uD[k]).mid())) for k in range(n))
print("max |J^T Z − Du| =", r)
zt = 3.0
et = acb(zt, 0)
s1 = list(acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs()); s1 += [acb(0)]*(n-len(s1))
kz = [et ** (-(j + 1)) for j in range(n)]
s2 = list((mser * acb_poly(kz)).truncate(n).coeffs()); s2 += [acb(0)]*(n-len(s2))
Jk = [sum((J[j,i]*kz[i] for i in range(n)), acb(0)) for j in range(n)]
g2 = [s1[j]+s2[j]-Jk[j] for j in range(n)]
print("max|g2_j| (unskaliert) =", max(float(abs(g).mid()) for g in g2), " g2_0 =", g2[0], " g2_5 =", g2[5])
# direkte Schwanzsumme: Σ_{k>N} ζ^{-k-1} [f^k]_{≤N} bis K
Pk = acb_poly([1])
for _ in range(n): Pk = (Pk * fpol).truncate(n)
tail = [acb(0)]*n
for k in range(n, 6*n):
    c = Pk.coeffs(); c = list(c) + [acb(0)]*(n-len(c))
    w = et ** (-(k+1))
    for j in range(n): tail[j] += w * c[j]
    Pk = (Pk * fpol).truncate(n)
print("max|tail_j| =", max(float(abs(g).mid()) for g in tail), " tail_0 =", tail[0], " tail_5 =", tail[5])
