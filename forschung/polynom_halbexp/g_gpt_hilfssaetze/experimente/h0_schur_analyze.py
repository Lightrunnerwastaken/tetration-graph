"""Summarize all-degree scans; plots and JSON have numerical status only."""
import json,math,argparse,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def describe(cases):
    n=np.array([c['n'] for c in cases]);sqrt=np.sqrt(n)
    result={'range':[int(n[0]),int(n[-1])],'count':len(n)}
    segments=[]
    for c in cases:
        if not segments or segments[-1]['point']!=c['fixed_argmax']:
            segments.append({'point':c['fixed_argmax'],'from_n':c['n'],'to_n':c['n']})
        else:segments[-1]['to_n']=c['n']
    result['fixed_argmax_segments']=segments
    for key in ['fixed_max','circle_sample_max','circle_taylor_majorant']:
        v=np.array([c[key] for c in cases]);s=sqrt*v
        indices=np.where((s[1:-1]>s[:-2])&(s[1:-1]>s[2:]))[0]+1
        result[key]={'max':float(max(v)),'max_n':int(n[np.argmax(v)]),
                     'scaled_max':float(max(s)),'scaled_max_n':int(n[np.argmax(s)]),
                     'scaled_min':float(min(s)),'scaled_min_n':int(n[np.argmin(s)]),
                     'scaled_last':float(s[-1]),'local_maxima':[[int(n[i]),float(s[i])] for i in indices],
                     'largest_up_step_relative':float(max(s[1:]/s[:-1]-1)),
                     'largest_down_step_relative':float(min(s[1:]/s[:-1]-1)),
                     'monotone_decreasing_unscaled':bool(np.all(np.diff(v)<0)),
                     'monotone_decreasing_scaled':bool(np.all(np.diff(s)<0))}
    for key in ['h_ratio','Q_norm','relative_pivot_2norm','h_over_diagonal']:
        v=np.array([c[key] for c in cases])
        result[key]={'min':float(min(v)),'min_n':int(n[np.argmin(v)]),
                     'max':float(max(v)),'max_n':int(n[np.argmax(v)]),'last':float(v[-1])}
    h=np.array([float(c['h']) for c in cases])
    result['h']={'min':float(min(h)),'max':float(max(h)),'positive_everywhere':bool(np.all(h>0)),
                 'strictly_increasing':bool(np.all(np.diff(h)>0))}
    q=np.array([c['Q_norm'] for c in cases])
    result['Q_growth_ratio']={'min':float(min(q[1:]/q[:-1])),'max':float(max(q[1:]/q[:-1]))}
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dir',type=Path,required=True)
    a=p.parse_args();records={}
    for name in ['kubisch','kneser_z']:
        path=a.dir/f'h0_schur_{name}.json'
        if path.exists():records[name]=json.loads(path.read_text(encoding='utf-8'))
    output={'status':'numerical finite-degree summary, not an asymptotic or interval statement','models':{}}
    fig,axs=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    colors={'kubisch':'#145a8d','kneser_z':'#b45824'}
    for name,d in records.items():
        cs=d['cases'];n=np.array([c['n'] for c in cs]);sq=np.sqrt(n);color=colors[name]
        summary={'completed':d['completed'],'all':describe(cs),'from_200':describe([c for c in cs if c['n']>=200]),
                 'sessions_seconds':d.get('sessions_seconds',[])+[d['seconds']]}
        directpath=a.dir/f'h0_schur_direct_{name}.json'
        if directpath.exists():
            direct=json.loads(directpath.read_text(encoding='utf-8'))
            summary['direct_plus300']={'degrees':[c['n'] for c in direct['cases']],
                                      'max_point_difference':max(c['max_increment_difference_scaled'] for c in direct['cases']),
                                      'max_circle_majorant_difference':max(c['majorant_difference_scaled'] for c in direct['cases']),
                                      'sessions':direct.get('sessions',[])}
        checkpath=a.dir/f'h0_schur_{name}_plus300.json'
        if checkpath.exists():
            dd=json.loads(checkpath.read_text(encoding='utf-8'));by={c['n']:c for c in dd['cases']}
            diff=[abs(v-by[c['n']]['points'][k])/(1+abs(v)) for c in cs if c['n'] in by for k,v in c['points'].items()]
            summary['plus300']={'last_n':dd['last_n'],'completed':dd['completed'],
                                'max_difference_published_points':max(diff),
                                'max_Q_relative_difference':max(abs(c['Q_norm']/by[c['n']]['Q_norm']-1) for c in cs if c['n'] in by),
                                'max_h_relative_difference_float':max(abs(float(c['h'])/float(by[c['n']]['h'])-1) for c in cs if c['n'] in by)}
        output['models'][name]=summary
        label='kubisch' if name=='kubisch' else 'Kneser'
        axs[0,0].plot(n,sq*np.array([c['fixed_max'] for c in cs]),color=color,label=label+': γ*-Proben')
        axs[0,0].plot(n,sq*np.array([c['circle_sample_max'] for c in cs]),color=color,ls='--',label=label+': Kreisraster')
        axs[0,1].plot(n,sq*np.array([c['circle_taylor_majorant'] for c in cs]),color=color,label=label)
        axs[1,0].plot(n,[c['h_ratio'] for c in cs],color=color,label=label)
        axs[1,1].plot(n,np.log10([float(c['h']) for c in cs]),color=color,label=label)
        axs[2,0].plot(n,np.log10([c['Q_norm'] for c in cs]),color=color,label=label)
        axs[2,1].plot(n,np.log10([c['relative_pivot_2norm'] for c in cs]),color=color,label=label)
    axs[0,0].set_title('√n · maximaler ausgewerteter Zuwachs (numerisch)')
    axs[0,1].set_title('√n · volle Kreis-Taylormajorante (numerisch)')
    axs[0,1].axhline(3.13,color='gray',ls=':',label='3,13: kubischer Kandidat ab n=200')
    axs[1,0].set_title('Schur-Pivotquotient hₙ / hₙ₋₁')
    axs[1,0].yaxis.set_major_formatter(FormatStrFormatter('%.6f'))
    axs[1,1].set_title('log₁₀ hₙ (alle gemessenen hₙ > 0)')
    axs[2,0].set_title('log₁₀ ‖Qₙ‖₂')
    axs[2,1].set_title('log₁₀ relative Pivotdiagnose (konservative 2-Norm)')
    for ax in axs.flat:ax.set_xlabel('Grad n');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.suptitle('H₀ Runde 4 — jeder Grad, Mittelpunktarithmetik; kein All-n-Zertifikat',fontsize=13)
    fig.savefig(a.dir/'H0_RUNDE4_KURVEN.png',dpi=160)
    (a.dir/'h0_schur_summary.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'summary_written':str(a.dir/'h0_schur_summary.json'),
                      'ranges':{k:v['all']['range'] for k,v in output['models'].items()}},ensure_ascii=False))

if __name__=='__main__':main()
