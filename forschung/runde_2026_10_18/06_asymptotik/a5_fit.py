import json, mpmath as mp
mp.mp.dps=40
g=[mp.mpf(s) for s in json.load(open("daten/g_minus2.json"))]
ns=list(range(40,141,10))
# Modell: (-1)^n g_n n log n = 1 + a/l + b/l^2 + c/l^3, l = log n  (Theorie: a = g0 - Euler-gamma)
A=mp.matrix([[1/mp.log(n),1/mp.log(n)**2,1/mp.log(n)**3] for n in ns]); y=mp.matrix([(-1)**n*g[n]*n*mp.log(n)-1 for n in ns])
sol=mp.lu_solve(A.T*A,A.T*y)
print("Fit a,b,c:",[mp.nstr(s,6) for s in sol]," Theorie a = g0-gamma =",mp.nstr(g[0]-mp.euler,6))
