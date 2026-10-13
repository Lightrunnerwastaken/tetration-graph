"""Analytic normalized summation and its fractional kernel (numerical audit).

No height polynomial interpolation or shared tetration runtime is used.
"""
import json,time,argparse
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).absolute().parent
def B(nu,z):
    if z==0:return mp.mpf(0)
    return mp.gamma(z+nu)*mp.rgamma(z)/mp.gamma(nu+1)
def ap_sum(f,z,T=None):
    f0=f(0)
    line=z*mp.quad(lambda u:f(z*u),[0,1])
    def integrand(t):
        if not t:return (mp.diff(f,0)-mp.diff(f,z))/mp.pi
        return 1j*(f(1j*t)-f(-1j*t)-f(z+1j*t)+f(z-1j*t))/mp.expm1(2*mp.pi*t)
    return line+(f(z)-f0)/2+mp.quad(integrand,[0,mp.mpf('.25'),1,3,8,mp.inf if T is None else T])
def kernel(s,u):return mp.gamma(u+s)*mp.rgamma(u+1)/mp.gamma(s)
def ap_fractional(f,s,z):
    f0=f(0)
    line=z*mp.quad(lambda u:kernel(s,z-z*u)*f(z*u),[0,1])
    def numerator(t):
        return (kernel(s,z-1j*t)*f(1j*t)-kernel(s,z+1j*t)*f(-1j*t)
            -kernel(s,-1j*t)*f(z+1j*t)+kernel(s,1j*t)*f(z-1j*t))
    def integrand(t):
        if not t:return 1j*mp.diff(numerator,0)/(2*mp.pi)
        return 1j*numerator(t)/mp.expm1(2*mp.pi*t)
    return line+(f(z)-kernel(s,z)*f0)/2+mp.quad(integrand,[0,mp.mpf('.25'),1,3,8,mp.inf])
def run(dps):
    mp.mp.dps=dps;start=time.monotonic();records=[]
    tests=[]
    v=lambda z:z*(z-1)
    w=lambda z:z*(z-1)*(z+1)/3
    q1=lambda z:z*(z-1)*(z-2)*(z+1)*(8*z+1)/60
    for z in [mp.mpf('.5'),mp.mpc('1.5','.5')]:
        tests += [('S_polynomial',z,lambda z=z:ap_sum(v,z),w(z)),
            ('Q1_polynomial',z,lambda z=z:ap_sum(lambda t:(w(t)-v(t))*(2*t-1),z),q1(z))]
        for nu in [mp.mpf('1.5'),mp.mpf('2.5')]:
            tests += [('S_gamma_'+str(nu),z,lambda z=z,nu=nu:ap_sum(lambda t:B(nu,t),z),B(nu+1,z)),
                ('Lhalf_gamma_'+str(nu),z,lambda z=z,nu=nu:ap_fractional(lambda t:B(nu,t),mp.mpf('.5'),z),B(nu+mp.mpf('.5'),z))]
        tests += [('Lhalf_v',z,lambda z=z:ap_fractional(v,mp.mpf('.5'),z),2*(B(mp.mpf('2.5'),z)-B(mp.mpf('1.5'),z)))]
    for name,z,fn,expected in tests:
        val=fn();err=abs(val-expected)
        records.append({'test':name,'z_re':str(z.real),'z_im':str(z.imag),
            'value_re':mp.nstr(val.real,dps-8),'value_im':mp.nstr(val.imag,dps-8),
            'expected_re':mp.nstr(expected.real,dps-8),'expected_im':mp.nstr(expected.imag,dps-8),
            'abs_error':mp.nstr(err,15)})
        print(name,z,mp.nstr(err,8),round(time.monotonic()-start,2),flush=True)
        if time.monotonic()-start>1080:raise TimeoutError('18 minute budget')
    result={'dps':dps,'seconds':time.monotonic()-start,'records':records}
    (HERE/f'data/abel_plana_{dps}.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dps',type=int,default=50);a=p.parse_args();run(a.dps)
