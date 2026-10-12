"""Summarize actual selected critical profiles separately from Gamma identities."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

HERE=Path(__file__).absolute().parent
mp.mp.dps=70
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,"axes.spines.right":False})

def value(item):
    z=item["value"]
    return mp.mpc(*[mp.mpf(k) for k in z]) if isinstance(z,list) else mp.mpc(mp.mpf(z))

profiles=[]
for r in range(4,9):
    choices=[json.loads(p.read_text(encoding="utf-8")) for p in HERE.glob(f"critical_profile_r{r}_*.json")]
    choices=[z for z in choices if z["status"]=="completed"]
    if choices:
        profiles.append(max(choices,key=lambda z:z["requested_dps"]))

rows=[];anchor_errors=[];fold_errors=[];recursion_errors=[]
for p in profiles:
    n=p["rank"];m=n-2;eta=mp.mpf(p["base"])
    for item in p["values"]:
        x=mp.mpf(item["x"]); y=value(item)
        if x <= 0 and x==int(x): anchor_errors.append(abs(y-x-1))
        if x==mp.mpf(p["fold_height"]):fold_errors.append(abs(y-x))
        if x not in [mp.mpf("-.5"),mp.mpf("-.25"),mp.mpf(".25"),mp.mpf(".5"),mp.mpf(1),mp.mpf("1.5"),mp.mpf(2),mp.mpf(3)]:continue
        b=mp.gamma(x+m)*mp.rgamma(x)/mp.gamma(m+1)
        d=y.real-x-1;a=d/b;u=m**(x-1)*a
        row={"rank":n,"m":m,"eta":eta,"x":x,"H":y.real,"imaginary_part":y.imag,"D":d,"B":b,"A":a,"U":u}
        if x==2:row["base_sensitivity_lower_bound"]=(4-y.real)/(2-eta)
        rows.append(row)
    for z in p.get("recursion_checks",{}).get("values",[]):
        recursion_errors.append(abs(value({"value":z})))

fields=list(rows[0])+["base_sensitivity_lower_bound"]
with (HERE/"critical_profile_audit.csv").open("w",encoding="utf-8",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    writer.writerows({k:mp.nstr(v,38) if isinstance(v,mp.mpf) else v for k,v in row.items()} for row in rows)
controls=[]
low_path=HERE/"critical_profile_r7_20.json";high_path=HERE/"critical_profile_r7_30.json"
if low_path.exists() and high_path.exists():
    lo=json.loads(low_path.read_text());hi=json.loads(high_path.read_text())
    controls=[abs(value(x)-value(y))/(1+abs(value(y))) for x,y in zip(lo["values"],hi["values"])]
meta={
    "profiles": [{k:p[k] for k in ["rank","base","requested_dps","status","seconds"]} for p in profiles],
    "all_completed_run_seconds": {p.name:json.loads(p.read_text())["seconds"] for p in HERE.glob("critical_profile_r*_*.json")},
    "max_anchor_residual":mp.nstr(max(anchor_errors),35),
    "max_fold_value_residual":mp.nstr(max(fold_errors),35),
    "max_sampled_recursion_residual":mp.nstr(max(recursion_errors),35),
    "rank7_max_relative_difference_20_vs_30_requested_dps":mp.nstr(max(controls),35),
    "accuracy_note":"Numerical consistency only. GP requested digits and printed digits are not rigorous error bounds. No interval certification. Rank 8 uses experimental two-arc selection.",
}
(HERE/"profile_metadata.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")

fig,ax=plt.subplots(figsize=(8.0,4.8),layout="constrained")
for p in profiles:
    points=sorted((float(mp.mpf(z["x"])),float(value(z).real)) for z in p["values"] if mp.mpf(z["x"])>=mp.mpf("-.5"))
    ax.plot(*zip(*points),marker="o",markersize=3.5,label=f"H{p['rank']} bei η{p['rank']}")
ax.plot([-.5,3.8],[.5,4.8],color="#666666",ls="--",lw=1,label="x + 1")
ax.plot([-.5,3.8],[-.5,3.8],color="#aaaaaa",ls=":",lw=1,label="x")
ax.set(xlabel="Höhe x",ylabel="H_N(η_N, x)",title="Kritische Familie: neue Stichproben, keine Grenzextrapolation",xlim=(-.55,3.85),ylim=(.4,4.85))
ax.legend(fontsize=9,ncol=2);ax.grid(alpha=.16)
for ext in ["png","svg"]:fig.savefig(HERE/f"kritische_profile.{ext}",dpi=180)
plt.close(fig)

with (HERE/"linear_half_step.csv").open(encoding="utf-8",newline="") as f: toy=list(csv.DictReader(f))
fig,ax=plt.subplots(figsize=(8.0,4.5),layout="constrained")
for col,label,color in [("v","Eingangsabweichung v(x)","#287aa9"),("half_step","Halber linearer Rangschritt","#d37b29"),("full_step","Voller linearer Rangschritt","#36865d")]:
    ax.plot([float(z["x"]) for z in toy],[float(z[col]) for z in toy],label=label,lw=2,color=color)
ax.axhline(0,color="#777777",lw=.7);ax.axvline(1,color="#cccccc",lw=.7)
ax.set(xlabel="Höhe x",ylabel="Lineare Abweichung vom Successor",title="Exakter halber Schritt im linearen Gamma-Modell")
ax.legend(fontsize=9);ax.grid(alpha=.15)
for ext in ["png","svg"]:fig.savefig(HERE/f"linearer_halbschritt.{ext}",dpi=180)
plt.close(fig)

print(json.dumps(meta,indent=2))
for p in profiles:
    r=p["rank"];sample={str(mp.mpf(z["x"])):value(z).real for z in p["values"]}
    lower=(4-sample["2.0"])/(2-mp.mpf(p["base"]))
    print(r,"F(-.5)",mp.nstr(sample["-0.5"],12),"F(.5)",mp.nstr(sample["0.5"],12),"F(2)",mp.nstr(sample["2.0"],12),"K lower",mp.nstr(lower,10))
