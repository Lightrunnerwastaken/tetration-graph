import json, mpmath as mp
mp.mp.dps=80
B=[mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))]
A=[mp.mpf(s) for s in json.load(open("daten/taylor_0_1_512_60.json"))]
def d(c,n): return c[n]-(-1)**(n+1)/(n*mp.mpf(2)**n)
for n in range(1,150,4):
    print(n, mp.nstr(d(B,n)*(-3)**n,10), mp.nstr(d(B,n+1)*(-3)**(n+1),10), mp.nstr(d(A,n)*(-3)**n,10))
