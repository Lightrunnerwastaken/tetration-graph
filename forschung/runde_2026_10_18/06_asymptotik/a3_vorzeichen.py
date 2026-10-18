import json, mpmath as mp
mp.mp.dps=60
B=[mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))]
A=[mp.mpf(s) for s in json.load(open("daten/taylor_0_1_512_60.json"))]
s=''.join('+' if B[n]>0 else '-' for n in range(0,150))
print(s)
neg=[n for n in range(150) if B[n]<0]; print("erste negative:",neg[:5])
alt=[n for n in range(1,150) if mp.sign(B[n])!=(-1)**(n+1)]
print("letzte Abweichung vom Muster (-1)^(n+1):", alt[-5:])
for n in range(20,30): print(n, mp.nstr(A[n],20), mp.nstr(B[n],20))
