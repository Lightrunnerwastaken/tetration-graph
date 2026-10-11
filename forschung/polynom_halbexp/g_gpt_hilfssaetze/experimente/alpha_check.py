"""Interval check of alpha < 0.37, conditional on the upstream certificate."""
import argparse
import json
import sys
from pathlib import Path
from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
PACKAGE = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PACKAGE.parent)
    parser.add_argument("--output-dir", type=Path, default=PACKAGE / "ergebnisse")
    args = parser.parse_args()
    data = json.loads((args.project_root / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))
    ctx.prec = 512
    midpoint, error = data["a"][0][:2]
    # Enlarge the upstream error to include rounding of its decimal strings.
    a0 = arb(midpoint, arb(error) + arb("1e-139"))
    rho = abs(-acb(-1).lambertw(-1))
    alpha = a0 / rho
    assert alpha.is_finite() and alpha > 0 and alpha < arb("0.37")
    result = {
        "status": "rigorous interval check conditional on the upstream coefficient certificate; midpoint rounding enlarged",
        "alpha_interval": alpha.str(40),
        "upper_bound_037_verified": True,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "alpha_check.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
