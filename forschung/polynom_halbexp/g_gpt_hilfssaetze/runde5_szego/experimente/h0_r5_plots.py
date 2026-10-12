"""Publication-ready finite numerical diagnostics for Runde 5."""
import json,math
from pathlib import Path
from h0_r5_paths import locations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
folder,ROUND4,ROOT=locations(__file__)
families=json.loads((folder/'family_analysis.json').read_text())['families']
old=json.loads((folder/'old_asymptotics.json').read_text())
profiles=json.loads((folder/'profile_analysis.json').read_text())
spiral=json.loads((folder/'spiral_analysis.json').read_text())
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
fig,axs=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
def lam(d):return float(next(v for v in d['fits'] if v['low']==250 and v['powers']==[0,2,3,4])['lambda'])
jet=sorted([(float(d['name'][4:]),lam(d)) for d in families if d['name'].startswith('jet_')]+[(0,lam(next(d for d in families if d['name']=='base')))])
axs[0,0].plot(*zip(*jet),'o-');axs[0,0].set(xlabel='ε bei festen L und μ',ylabel='λ, extrapoliert',title='Globale Daten ändern die Hauptordnung')
axs[0,0].set_xticks([-.01,-.005,0,.005,.01])
tau=sorted([(float(d['name'][4:]),lam(d)) for d in families if d['name'].startswith('tau_')]+[(1,lam(next(d for d in families if d['name']=='base')))])
base_lambda=lam(next(d for d in families if d['name']=='base'))
axs[0,1].semilogx([v[0] for v in tau],[(v[1]-base_lambda)*1e10 for v in tau],'o-');axs[0,1].set(xlabel='M-Amplitude τ > 0',ylabel='(λ_fit − λ_fit,τ=1) × 10¹⁰',title='Positive M-Amplituden: gleiche Hauptordnung?')
for model,label in [('kubisch','kubisch'),('kneser_z','Kneser')]:
    d=json.loads((ROUND4/f'h0_schur_{model}.json').read_text());cases=d['cases'];ns=np.array([c['n'] for c in cases[1:]])
    ratios=np.array([math.log(float(cases[j+1]['h'])/float(cases[j]['h'])) for j in range(len(cases)-1)])
    low=800 if model=='kubisch' else 450
    f=next(v for v in old[model]['ratio_fits'] if v['low']==low and v['powers']==[0,2,3,4])
    axs[0,2].plot(ns,ns**2*(ratios-float(f['log_coefficients'][0])),label=label)
    axs[1,0].plot([c['n'] for c in cases],[math.log(c['Q_norm'])-math.log(float(c['h'])) for c in cases],label=label)
    a=profiles[model]['corner_scaled_absolute'];axs[1,2].plot([int(n) for n in a],[a[n] for n in a],'.-',label=label)
axs[0,2].set(xlabel='n',ylabel='n² (log rₙ − log λ_fit)',title='Numerischer Term der Ordnung 1/n²',ylim=(-.48,-.41));axs[0,2].legend()
axs[1,0].set(xlabel='n',ylabel='log(‖Qₙ‖/hₙ)',title='Exponentieller Normquotient');axs[1,0].legend()
diffs=spiral['models']['kubisch']['differences'];diffs=sorted([v for v in diffs if v['m']==2*v['n']],key=lambda v:v['n'])
axs[1,1].loglog([v['n'] for v in diffs],[v['fixed_angle_normalized_W_difference'] for v in diffs],'o-',label='feste Winkel')
axs[1,1].loglog([v['n'] for v in diffs],[v['natural_normalized_W_difference'] for v in diffs],'o-',label='komplexe Kegelskala')
axs[1,1].set(xlabel='n, Vergleich n → 2n',ylabel='Wiener-Abstand der normierten Profile',title='Drehung entfernen verbessert die Annäherung');axs[1,1].legend()
axs[1,2].set(xlabel='n',ylabel='√n |Dₙ(L)|',title='Amplitude: Grenzwert weiterhin unsicher');axs[1,2].legend()
fig.suptitle('Runde 5: endliche Numerik; keine All-n-Schranken',fontsize=15)
fig.savefig(folder/'H0_RUNDE5_ASYMPTOTIK.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(12,10),constrained_layout=True)
angles=np.linspace(0,2*np.pi,501);z=np.exp(1j*angles)
by={v['n']:v for v in spiral['models']['kubisch']['profiles']}
for n in [200,400,800,1000]:
    c=np.array([complex(*map(float,v)) for v in by[n]['natural_shape_coefficients_first_40']]);p=np.polynomial.polynomial.polyval(z,c)
    axs[0,0].plot(p.real,p.imag,label=f'n={n}')
axs[0,0].set(xlabel='Re Sₙ(ζ)',ylabel='Im Sₙ(ζ)',title='Normiertes Profil auf der komplexen Kegelskala');axs[0,0].legend();axs[0,0].axis('equal')
phase_data=[]
for file in folder.glob('profiles_kubisch_*.json'):
    for c in json.loads(file.read_text())['cases']:
        if c['n'] in [200,300,400,500,600,700,800,900,1000]:
            aa=np.array([complex(*map(float,v)) for v in c['w_circle_taylor_scaled'][:2]]);phase_data.append((c['n'],np.angle(aa[1]/aa[0])))
phase_data=sorted(set(phase_data));nn=np.array([v[0] for v in phase_data]);ph=np.array([v[1] for v in phase_data]);sigma=float(spiral['inverse_kappa'][1])
axs[0,1].plot(np.log(nn),ph,'o-',label='Phase des ersten Taylorquotienten')
axs[0,1].plot(np.log(nn),ph[-1]+sigma*np.log(nn/nn[-1]),'--',label='Steigung Im(1/κ)')
axs[0,1].set(xlabel='log n',ylabel='Phase in Radiant',title='Langsame räumliche Drehung');axs[0,1].legend()
for v in diffs:
    axs[1,0].plot([v['n']],[v['natural_normalized_W_difference']],'o',color='C0')
axs[1,0].plot([v['n'] for v in diffs],[v['natural_normalized_W_difference'] for v in diffs],'-',color='C0')
axs[1,0].set(xlabel='n',ylabel='‖S₂ₙ − Sₙ‖W',title='Numerische Annäherung nach Normierung')
pp=json.loads((folder/'phase_analysis.json').read_text())['points'];labels=list(pp);omega=[float(pp[k]['omega']) for k in labels]
axs[1,1].bar(labels,omega);axs[1,1].set(ylabel='ω(x), geschätzter Träger',title='Verschiedene gemessene Phasenträger auf γ*',ylim=(-1.55,-1.42))
fig.suptitle('Komplexe Profile: endlich ausgewertet, Grenzprofil nur Kandidat',fontsize=15)
fig.savefig(folder/'H0_RUNDE5_PROFILE.png',dpi=160);plt.close(fig)
print('Two plots saved.')
