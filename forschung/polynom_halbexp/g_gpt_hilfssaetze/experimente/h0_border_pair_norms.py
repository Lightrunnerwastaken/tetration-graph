"""Read existing W inverses: norm of the rank-one degree increment operator.

Float64 diagnostic using already-computed inverses and singular values;
no matrix solve and no fresh SVD. Not a certificate or asymptotic theorem.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    result={'status':'float64 diagnostic from existing w01 inverses; no new SVD or interval certificate','cases':[]}
    for N in [200,300,400,500,600,700,800,900,1000]:
        z=np.load(a.project_root/f'w_grenzschicht/w01_N{N}.npz')
        G=z['G']
        qnorm=np.linalg.norm(G[:,-1])/abs(G[-1,-1])
        dn_coeff=qnorm*np.linalg.norm(G[-1,:])
        sig=float(z['S'][0])
        result['cases'].append({'N':N,'Q_norm':float(qnorm),
                                'degree_increment_operator_norm':float(dn_coeff),
                                'existing_W_inverse_sigma1':sig,
                                'ratio_to_existing_sigma1':float(dn_coeff/sig)})
    ns=np.array([x['N'] for x in result['cases']])
    ln=np.log([x['degree_increment_operator_norm'] for x in result['cases']])
    fit=np.linalg.lstsq(np.column_stack([ns,np.log(ns),np.ones(len(ns))]),ln,rcond=None)[0]
    result['preasymptotic_fit']={'linear_N':float(fit[0]),'log_N':float(fit[1]),'constant':float(fit[2]),
                                'warning':'finite-window fit, not a determination of the asymptotic rate c'}
    result['later_window_fits']=[]
    for lower in [400,600]:
        mask=ns>=lower
        f2=np.linalg.lstsq(np.column_stack([ns[mask],np.log(ns[mask]),np.ones(np.count_nonzero(mask))]),ln[mask],rcond=None)[0]
        result['later_window_fits'].append({'lower_N':lower,'upper_N':1000,'linear_N':float(f2[0]),'log_N':float(f2[1]),'constant':float(f2[2])})
    a.output_dir.mkdir(parents=True,exist_ok=True)
    (a.output_dir/'h0_border_pair_norms.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
