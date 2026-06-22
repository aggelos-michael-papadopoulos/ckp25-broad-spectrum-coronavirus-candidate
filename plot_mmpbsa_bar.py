#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SYSTEMS = [
    ("MERS_B_Spike_md", "MERS-CoV\n(4KR0)", "#9B59B6"),
    ("hCoV_NL63_E_Spike_md", "HCoV-NL63\n(3KBH)", "#E67E22"),
    ("SARS-CoV-1_E_spike_md", "SARS-CoV-1\n(2AJF)", "#3498DB"),
    ("hCoV_HKU1_A_spike_md", "HCoV-HKU1\n(8Y7Y)", "#1ABC9C"),
    ("hCoV_229E_E_spike_md", "HCoV-229E\n(6ATK)", "#2ECC71"),
]


def parse_delta_total(path: Path) -> tuple[float, float]:
    for line in path.read_text().splitlines():
        if line.strip().startswith("DELTA TOTAL"):
            parts = line.split()
            return float(parts[2]), float(parts[3])
    raise ValueError(f"DELTA TOTAL not found in {path}")


def main():
    parser = argparse.ArgumentParser(description="Plot MM-GBSA binding free energies.")
    parser.add_argument("--results-dir", default="runs")
    parser.add_argument("--out-dir", default="paper_figures")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = []
    delta_g = []
    std_dev = []
    colors = []

    for folder, label, color in SYSTEMS:
        result_file = results_dir / folder / "FINAL_RESULTS_MMPBSA.dat"
        if not result_file.exists():
            print(f"Skipping {folder}: {result_file} not found.")
            continue
        mean, sd = parse_delta_total(result_file)
        labels.append(label)
        delta_g.append(mean)
        std_dev.append(sd)
        colors.append(color)

    if not delta_g:
        raise SystemExit(f"No FINAL_RESULTS_MMPBSA.dat files found in {results_dir}")

    plt.style.use("default")
    plt.rcParams.update({
        "font.size": 12,
        "axes.linewidth": 1.2,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
    })

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axhline(0, color="black", linewidth=1.0, zorder=1)

    bars = ax.bar(
        labels,
        delta_g,
        yerr=std_dev,
        color=colors,
        edgecolor="none",
        width=0.5,
        capsize=6,
        error_kw=dict(ecolor="#2C3E50", lw=1.5, capthick=1.5),
        zorder=3,
    )

    for bar, val in zip(bars, delta_g):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            3.0,
            f"{val:.1f}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color="#2C3E50",
        )

    ax.set_ylabel(r"$\Delta$G$_{bind}$ (kcal/mol)", fontsize=13, fontweight="bold")
    ax.set_title(r"MM-GBSA binding free energy ($\Delta$G$_{bind}$)", fontsize=14, fontweight="bold", pad=15)
    ax.set_ylim(min(delta_g) - 12, 10)
    ax.set_yticks(np.arange(-50, 1, 10))
    ax.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = out_dir / "md_gbsa_bar.jpg"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved MM-GBSA bar chart to {out_path}")


if __name__ == "__main__":
    main()
