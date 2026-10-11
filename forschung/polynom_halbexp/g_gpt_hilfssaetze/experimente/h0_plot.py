"""Plot the saved H0 diagnostics without recomputing experiments."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    base = args.input_dir or package / "ergebnisse"
    output = args.output or package / "ergebnisse" / "H0_PRUEFBILD.png"
    cubic = json.loads((base / "h0_cauchy_kubisch.json").read_text(encoding="utf-8"))["cases"]
    kneser = json.loads((base / "h0_cauchy_kneser_z.json").read_text(encoding="utf-8"))["cases"]
    integral = json.loads((base / "h0_kernel_probe.json").read_text(encoding="utf-8"))["cases"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.9), layout="constrained")
    for data, label, color, marker in [(cubic,"kubisch", "#235e91", "o"), (kneser,"Kneser", "#c46e2d", "x")]:
        ns = [case["N"] for case in data]
        axes[0].plot(ns, [case["G_whole_circle_over_sqrtN"] for case in data], marker=marker, color=color, label=f"G: Eckkreis, {label}")
        axes[0].plot(ns, [case["gamma_sample_G_over_sqrtN"] for case in data], marker=marker, color=color, linestyle="--", label=f"G: feste γ*-Proben, {label}")
    axes[0].set(xlabel="N", ylabel="Majorante / √N", title="SVD-freie Cauchy-Antworten", ylim=(0.35, 0.9))
    axes[0].legend(fontsize=8.5)
    ns = [case["N"] for case in integral]
    for field, label, color, marker in [("E_abs", "Gemeinsames |I_f − B_y|", "#2b7550", "o"), ("absolute_integral", "Integral der Beträge", "#b54a47", "s"), ("weighted_CS_bound", "Gewichtete CS-Majorante", "#8a62a0", "^")]:
        axes[1].semilogy(ns, [case[field] for case in integral], color=color, marker=marker, label=label)
    axes[1].set(xlabel="N", ylabel="Betrag", title="Auslöschung im gemeinsamen Integral")
    axes[1].legend(fontsize=8.5)
    for ax in axes:
        ax.grid(True, alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Numerische Diagnostik von (H₀) · 11.10.2026", fontsize=14)
    fig.supxlabel("Stichpunkte, keine All-N-Schranke. Rechts: approximative W-Eckvektoren, Auswertung mit 2400 Bit.", fontsize=8.5)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    print(output)

if __name__ == "__main__":
    main()
