"""Tables and standalone figures for the analytic half-rank round."""
import json,csv,math,time
from pathlib import Path
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).absolute().parent
def read(name):return json.loads((HERE/'data'/name).read_text(encoding='utf-8'))
def run():
    mp.mp.dps=100
    raw=read('direct_coefficients_K512_500.json');a=np.array([float(mp.mpf(r['value'])) for r in raw['coefficients']])
    samples=read('evaluation_K512_500.json')['samples']
    rank=read('rank_newton_N128_400.json')
    rankrows=[]
    for point in ['1/2','3/2','5/2','3']:
        from fractions import Fraction as F
        data=[F(r['values'][point]) for r in rank['integer_rank']];co=[]
        while data:
            co.append(data[0]);data=[d-c for c,d in zip(data,data[1:])]
        total=mp.mpf(0);bj=mp.mpf(1)
        for j,c in enumerate(co):
            if j:bj*= (mp.mpf('.5')-j+1)/j
            total+=bj*mp.mpf(c.numerator)/c.denominator
            if j in [16,32,64,96,128]:rankrows.append({'x':point,'cutoff':j,'q':mp.nstr(total,60)})
    with (HERE/'data/rank_first_evaluation.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['x','cutoff','q']);w.writeheader();w.writerows(rankrows)
    with (HERE/'data/evaluation.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(samples[0]));w.writeheader();w.writerows(samples)
    old=list(csv.DictReader((HERE/'inputs/model_size_audit.csv').open(encoding='utf-8')))
    new=[r for r in samples if r['x_re']=='0.5' and r['x_im']=='0.0']
    plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':.25,'figure.dpi':150})
    fig,axs=plt.subplots(1,2,figsize=(12,4.1))
    axs[0].plot([int(r['K']) for r in old],[float(r['Q_half_x0p5']) for r in old],'o-',label='alte endliche Polynomrekonstruktion')
    axs[0].plot([r['cutoff'] for r in new],[float(r['q_re']) for r in new],'o-',label='intrinsische Höhen-Newtonreihe')
    rr=[r for r in rankrows if r['x']=='1/2']
    axs[0].plot([r['cutoff'] for r in rr],[float(r['q']) for r in rr],'s--',label='unabhängig: zuerst im Rang')
    axs[0].set_xscale('log',base=2);axs[0].set_xlabel('K bzw. Anzahl der Terme');axs[0].set_ylabel('quadratischer Halbrangwert bei x=0.5')
    axs[0].set_ylim(-.055,.065);axs[0].legend(fontsize=8)
    j=np.arange(8,513);scaled=(-1.)**j*a[j]*j**1.5
    A=float(-1/(2*mp.pi*mp.sqrt(mp.pi)))
    H=mp.digamma(mp.mpf('3.5'))+mp.euler;beta=(mp.digamma(mp.mpf('2.25'))-mp.digamma(mp.mpf('1.75')))/2
    B=float((-2*mp.digamma(-mp.mpf('.5'))/mp.pi-2*(H+beta)/mp.pi)/mp.gamma(-mp.mpf('.5')))
    axs[1].plot(j,scaled,label='gemessene Koeffizienten')
    # Only the additive constant is fitted; leading log terms are derived candidates.
    C=float(np.mean(scaled[-128:]-A*np.log(j[-128:])**2-B*np.log(j[-128:])))
    axs[1].plot(j,A*np.log(j)**2+B*np.log(j)+C,'--',label='vermutete Log-Asymptotik, C angepasst')
    axs[1].set_xlabel('j');axs[1].set_ylabel(r'$(-1)^j j^{3/2}\alpha_j$');axs[1].legend(fontsize=8)
    fig.suptitle('Quadratischer Halbrangjet: exakte Koeffizienten, numerische Partialsummen')
    fig.tight_layout();fig.savefig(HERE/'figures/konvergenz.png');fig.savefig(HERE/'figures/konvergenz.pdf');plt.close(fig)
    x=np.linspace(0,3,601);total=np.zeros_like(x);bj=np.ones_like(x);profiles={}
    for j,c in enumerate(a):
        if j:bj*= (x-j+1)/j
        total+=c*bj
        if j in [64,128,256,512]:profiles[j]=total.copy()
    fig,ax=plt.subplots(figsize=(8,4))
    for k,vals in profiles.items():ax.plot(x,vals,label=f'{k} Terme',lw=1.5)
    ax.set_xlabel('Höhe x');ax.set_ylabel(r'$q_{1/2}(x)$');ax.set_title('Analytischer quadratischer Testjet für v=x(x−1)')
    ax.legend();fig.tight_layout();fig.savefig(HERE/'figures/profil.png');fig.savefig(HERE/'figures/profil.pdf');plt.close(fig)
    # Preserve a machine-readable proof/numerics distinction.
    summary={'theorem_domain':'Re x > -1','coefficient_bound':'1500*(1+log(j))/j for j>=3',
        'partial_sum_at_half_512':next(r['q_re'] for r in new if r['cutoff']==512),
        'asymptotic_status':'sketch: uniform complex endpoint expansion and transfer remainder still open',
        'asymptotic_log2_coefficient':A,'asymptotic_log_coefficient':B,'fitted_constant_last128':C,
        'full_H3p5_constructed':False,'eta3p5_computed':False,
        'run_seconds':{p.name:json.loads(p.read_text(encoding='utf-8')).get('seconds') for p in (HERE/'data').glob('*.json') if 'seconds' in json.loads(p.read_text(encoding='utf-8'))}}
    (HERE/'data/analysis.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
