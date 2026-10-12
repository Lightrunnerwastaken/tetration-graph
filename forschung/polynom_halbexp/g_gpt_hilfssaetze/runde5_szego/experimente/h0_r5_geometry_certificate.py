"""Small scalar interval certificate for the spiral-circle geometry.

This does not certify any Schur inverse or limiting profile.
"""
import json,sys
from pathlib import Path
from h0_r5_paths import locations
from flint import arb,acb,ctx
sys.stdout.reconfigure(encoding='utf-8');ctx.prec=384
re='0.31813150520476413531265425158766451720351761387139986692237'
im='1.33723570143068940890116214319371061253950213846051241687'
c=acb(arb(re),arb(im));r=arb('1e-50')
lhs=abs(c.exp()-c)+(c.real+r).exp()*r*r/2
rhs=abs(c.exp()-1)*r
assert lhs<rhs
L=acb(arb(re)+arb('0 +/- 1e-50'),arb(im)+arb('0 +/- 1e-50'))
mu=L.sqrt();kappa=acb(0,2*arb.pi())/(acb(0,2*arb.pi())-mu.log());a=1/kappa
beta=1/arb('1.1182');diff=beta-a.real
ratio=((a.real-beta)*arb(50).log()).exp()/(1-(-beta*arb(50).log()).exp())
assert diff>0 and ratio<arb('1.05')
out={'status':'rigorous scalar interval certificate by Rouche and Arb; no matrix inverse certificate',
     'precision_bits':384,'root_disk_center':[re,im],'root_disk_radius':'1e-50',
     'Rouche_left_upper':str(lhs.upper()),'Rouche_right_lower':str(rhs.lower()),
     'mu_modulus':str(abs(mu)),'kappa':[str(kappa.real),str(kappa.imag)],
     'inverse_kappa':[str(a.real),str(a.imag)],'beta_minus_Re_inverse_kappa':str(diff),
     'circle_inclusion_ratio_n50':str(ratio),'fixed_enlargement_radius':'1.05',
     'uniform_degrees':'n>=50, by strict monotonicity of n^(Re(1/kappa)-beta)/(1-n^-beta)'}
folder,ROUND4,ROOT=locations(__file__)
(folder/'spiral_geometry_certificate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
