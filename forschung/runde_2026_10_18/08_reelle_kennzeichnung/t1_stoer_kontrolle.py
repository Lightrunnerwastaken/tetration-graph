import mpmath as mp, rk_lib as R, stoerungen as S
R.init(50); N=12
for name in ["sin:0.03","cos:0.01","moeb:0.05:0","moeb:0.2:0.25","gb:c"]:
    st=S.katalog(name); x=mp.mpf('0.137')
    ser=R.phi_reihe(st,x,N)
    f=lambda z: z+S.p_komplex(st,z)[0]
    # Taylor über Cauchy-Integral als Gegenprobe
    tay=mp.taylor(f,x,6,method='quad',radius=mp.mpf('0.1'))
    print(name, max(abs(ser[k]-tay[k]) for k in range(6)), mp.nstr(S.p_komplex(st,mp.mpf(0))[0],5))
