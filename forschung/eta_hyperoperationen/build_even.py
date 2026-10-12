"""Namespace a local copy of the existing two-arc Abel solver for rank 8.

No claim of a universal Kneser theorem: this creates an experimental,
selection-dependent extension whose residuals must be measured.
"""
import hashlib,json,re
import argparse
from common import HERE,snapshot
source=snapshot().read_text(encoding="utf-8")
p=argparse.ArgumentParser();p.add_argument("--pilot",action="store_true")
p.add_argument("--target",type=int,default=10);p.add_argument("--output",default="even8.gp");args=p.parse_args()
segment=source[source.index("knlsq(M, b) = {"):source.index("/* ecalle routine")]
symbols=set(re.findall(r"\bkn[A-Za-z0-9_]*\b",segment))|{"kp","krs","ksig","klam","klog","kper"}
replacement={name:"o8"+name for name in symbols}
replacement.update({"hexreg":"o8reg","hexerr":"o8err","pent":"hmap7","ipent":"hinv7"})
regex=re.compile(r"\b(?:"+"|".join(map(re.escape,sorted(replacement,key=len,reverse=True)))+r")\b")
text=regex.sub(lambda m:replacement[m.group()],segment)
# Limit target to the accuracy of the already constructed lower level.
text=text.replace("my(eps = 10.0^(-precis-8), target = precis - 3,",
                  "my(target = max(8,min(precis-10,floor(hrerror[7])-6)), eps = 10.0^(-target-6),")
if args.pilot:
    text=text.replace("target = max(8,min(precis-10,floor(hrerror[7])-6))",f"target = {args.target}")
    text=text.replace("getabstime() - tp > 400","getabstime() - tp > 10000")
    text=text.replace("alarm(1, o8knplan","alarm(10, o8knplan")
    text=text.replace("[0.45, 0.55, 0.65, 0.75, 0.85, 1.0, 1.15, 1.3]","[0.55, 0.65, 0.85, 1.0]")
    text=text.replace("for (ic = 1, 7,","for (ic = 1, 3,")
    text=text.replace("[0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5][ic]","[0, 0.75, 1.5][ic]")
prefix="""/* Experimental rank-8 copy of the upstream two-arc solver, separate state. */
o8reg=0;o8err=0;o8knforce=0;o8knlsqpass=2;o8knsolvep=2;
o8knreachtol=1.02;o8kndebug=1;o8knRmeas=0;o8knmaxk=250;
hmap7(z) = hval(7,z);
hinv7(z) = hinv(7,z);
"""
suffix="""
hcomplexinit(r)={
  if(r!=8,error("hcomplexinit: experimental solver currently supports rank 8 only"));
  if(!o8knhexinit(),error("hcomplexinit: no two-arc solution"));
  hrready[8]=1;hrkind[8]=2;hrp[8]=o8kp;hrlam[8]=o8klam;
  hrerror[8]=min(hrerror[7],o8err);hrtarget[8]=hrerror[8];
  return(hrerror[8]);
}
"""
(HERE/args.output).write_text(prefix+text+suffix,encoding="utf-8")
print(json.dumps({"cloned_symbols":len(replacement),"source_sha256":hashlib.sha256(snapshot().read_bytes()).hexdigest()}))
