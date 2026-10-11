"""Scientific plot of the recorded round-3 diagnostics; no new calculations."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(path): return json.loads(path.read_text(encoding='utf-8'))
def midpoint(text): return float(text.lstrip('[').split(' +/-')[0].rstrip(']'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    cube=load(a.data_dir/'h0_border_kubisch.json')['cases']
    pair=load(a.data_dir/'h0_border_pair_norms.json')['cases']
    cert=[]
    for name in ['h0_border_circle_cert.json','h0_border_circle_cert_remaining.json','h0_border_circle_cert_large.json','h0_border_circle_cert_1000.json']:
        if (a.data_dir/name).exists(): cert+=load(a.data_dir/name)['cases']
    cert=sorted(cert,key=lambda x:x['N'])
    assert [x['N'] for x in cert]==[200,300,400,600,800,1000]
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.7,4.5),layout='constrained')
    ns=[c['N'] for c in cube]
    axes[0].plot(ns,[c['fixed_gamma_max_sqrtN_increment'] for c in cube],'o-',label='feste γ*-Proben (numerisch)',color='#007d75')
    axes[0].plot(ns,[c['circle_max_sqrtN_increment'] for c in cube],'s-',label='24 Kreiswinkel (numerisch)',color='#365cb5')
    axes[0].plot([c['N'] for c in cert],[midpoint(c['sqrtN_whole_circle_increment_upper']) for c in cert],'+--',ms=10,lw=1,label='ganzer Kreis: Intervallmajorante',color='#a24a31')
    axes[0].axhline(3.13,color='#777777',ls=':',lw=1)
    axes[0].set(xlabel='Grad N',ylabel='√N · Norm des Grad-Zuwachses',title='Die ausgewerteten Zuwächse bleiben klein',ylim=(0,3.6))
    axes[0].legend(fontsize=9,loc='lower left')
    axes[0].grid(alpha=.15)
    axes[1].semilogy([c['N'] for c in pair],[c['degree_increment_operator_norm'] for c in pair],'o-',color='#7e4596',label='Grad-Zuwachs als Operator (W-Daten)')
    axes[1].semilogy([c['N'] for c in pair],[c['existing_W_inverse_sigma1'] for c in pair],'--',color='#555555',label='größter Singulärwert der Inversen (W-Daten)')
    axes[1].semilogy(ns,[c['fixed_gamma_max_sqrtN_increment']/np.sqrt(c['N']) for c in cube],'o-',color='#007d75',label='Grad-Zuwachs auf γ*-Proben')
    axes[1].set(xlabel='Grad N',ylabel='Norm / Betrag (logarithmische Skala)',title='Die Operatorverstärkung wächst im Messfenster')
    axes[1].legend(fontsize=9,loc='upper left')
    axes[1].grid(alpha=.15)
    fig.suptitle('Runde 3: exakte Grad-Erweiterung mit Quelle — keine All-N-Schranke',fontsize=14)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.output,dpi=180,facecolor='white')
    plt.close(fig)
    print(a.output)


if __name__=='__main__':
    main()
