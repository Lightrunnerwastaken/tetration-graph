"""Portable locations for local diagnostics and the published research folder."""
from pathlib import Path
import argparse

def locations(script):
    parent = Path(script).absolute().parents[1]
    if parent.name == "runde5_szego" and (parent / "ETAPPE1.md").exists():
        defaults=parent, parent.parent / "ergebnisse", parent.parent.parent
    else:
        defaults=(parent / "work/h0_round5", parent / "work/h0_round4",
                  Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp"))
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument("--data-dir",type=Path)
    parser.add_argument("--round4-dir",type=Path)
    parser.add_argument("--project-root",type=Path)
    args,_=parser.parse_known_args()
    overrides=(args.data_dir,args.round4_dir,args.project_root)
    return tuple(v.absolute() if v is not None else d for v,d in zip(overrides,defaults))
