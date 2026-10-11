"""Finite rigorous initial norm for the conditional degree-sum criterion."""
import argparse
import json
import sys
import time
from pathlib import Path
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
from h0_border_circle_cert import exact_model,inverse_circle_norm
from flint import arb, acb, acb_poly


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--N',type=int,default=199)
    p.add_argument('--radius-degree',type=int,default=200)
    a=p.parse_args()
    start=time.time()
    N=a.N
    bits=2*N+500
    J,L,rho,mu=exact_model(N,bits)
    Ji=J.inv() # checked full inverse
    beta=arb(1)/arb('1.1182')
    rad=rho*(-beta*arb(a.radius_degree).log()).exp()
    center=rho*L.conjugate()/(rho**2-rad**2)
    wrad=rho*rad/(rho**2-rad**2)
    upper,terms=inverse_circle_norm(Ji,center,wrad)
    scaled=upper/arb(N).sqrt()
    result={'status':'rigorous finite prefix norm on the whole corner disk via checked inverse and finite Taylor sum',
            'N':N,'radius_degree':a.radius_degree,'bits':bits,
            'C0_upper':str(upper.upper()),'C0_over_sqrtN_upper':str(scaled.upper()),
            'C0_over_sqrtN_less_than_6p26':bool(scaled<arb('6.26')),
            'Taylor_terms_first_12':[x.str(25) for x in terms[:12]],
            'Taylor_tail_after_12':sum(terms[12:],arb(0)).str(25),
            'seconds':time.time()-start}
    assert result['C0_over_sqrtN_less_than_6p26']
    a.output_dir.mkdir(parents=True,exist_ok=True)
    (a.output_dir/'h0_border_prefix_cert.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print('C0',upper.str(20),'C0/sqrtN',scaled.str(20),'seconds',result['seconds'],flush=True)


if __name__=='__main__':
    main()
