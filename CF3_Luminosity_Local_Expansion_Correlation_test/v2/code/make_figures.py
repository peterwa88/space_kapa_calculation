"""Recreate the manuscript figure and two supplementary sensitivity figures."""

from pathlib import Path
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cf3-v2")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 160,
    }
)


def binned(data: pd.DataFrame, nbins: int = 6, seed: int = 20260820) -> pd.DataFrame:
    """Return equal-count residual bins with bootstrap 95% intervals."""
    data = data.copy()
    data["bin"] = pd.qcut(data["resid_logL"], nbins, labels=False, duplicates="drop")
    rng = np.random.default_rng(seed)
    rows = []
    for _, subset in data.groupby("bin", observed=True):
        x = subset["resid_logL"].to_numpy()
        y = subset["resid_logH"].to_numpy()
        boot = np.empty(4000)
        for index in range(4000):
            ii = rng.integers(0, len(subset), len(subset))
            boot[index] = y[ii].mean()
        low, high = np.quantile(boot, [0.025, 0.975])
        rows.append(
            {
                "x": x.mean(),
                "y": 100.0 * (10.0 ** y.mean() - 1.0),
                "low": 100.0 * (10.0 ** low - 1.0),
                "high": 100.0 * (10.0 ** high - 1.0),
                "n": len(subset),
            }
        )
    return pd.DataFrame(rows)


def manuscript_figure() -> None:
    """Create the simplified primary-sample figure used in the final manuscript."""
    data = pd.read_csv(DATA / "analysis_primary.csv")
    summary = binned(data)

    # Frisch-Waugh-Lovell residual regression: this slope equals the adjusted
    # luminosity coefficient from the full nuisance model.
    beta, intercept = np.polyfit(data["resid_logL"], data["resid_logH"], 1)
    x_line = np.linspace(summary["x"].min() - 0.03, summary["x"].max() + 0.03, 250)
    y_line = 100.0 * (10.0 ** (intercept + beta * x_line) - 1.0)

    fig, ax = plt.subplots(figsize=(6.5, 3.75))
    ax.axhline(0, color="#7f8790", lw=0.9, ls="--", zorder=1)
    ax.plot(x_line, y_line, color="#1f5b99", lw=2.2,
            label="Adjusted positive trend", zorder=2)
    ax.errorbar(
        summary["x"], summary["y"],
        yerr=[summary["y"] - summary["low"], summary["high"] - summary["y"]],
        fmt="o", color="#173f5f", ecolor="#4b78a4", ms=6.0,
        capsize=3.5, capthick=1.0, lw=1.2, zorder=3,
        label="Six equal-count bins (95% bootstrap CI)",
    )
    ax.set_title("Positive luminosity-local-expansion-proxy association", fontweight="bold")
    ax.set_xlabel("Nuisance-adjusted log$_{10}$ group luminosity residual (dex)")
    ax.set_ylabel("Adjusted local-expansion-proxy residual (%)")
    ax.grid(axis="y", color="#d8dde3", lw=0.65, alpha=0.9)
    ax.text(
        0.03, 0.96,
        r"Primary sample: $N=858$; $\beta=0.0384$ dex/dex; $p=0.0176$",
        transform=ax.transAxes, va="top", color="#333333", fontsize=9.5,
    )
    ax.legend(frameon=False, loc="lower right", fontsize=8.8)
    fig.tight_layout()
    fig.savefig(OUT / "figure1_primary_positive_association.png", bbox_inches="tight", dpi=320)
    plt.close(fig)


def supplementary_binned_figure() -> None:
    """Compare the grouped-distance primary sample with the SN Ia-only sample."""
    samples = [
        ("Primary grouped-distance sample", "primary", "#1f5b99"),
        ("SN Ia-only distance sample", "sn_only", "#a65a2e"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharey=True)
    for ax, (title, name, color) in zip(axes, samples):
        data = pd.read_csv(DATA / f"analysis_{name}.csv")
        summary = binned(data, seed=20260819)
        ax.axhline(0, color="#7a7a7a", lw=0.9, ls="--", zorder=1)
        ax.errorbar(
            summary["x"], summary["y"],
            yerr=[summary["y"] - summary["low"], summary["high"] - summary["y"]],
            fmt="o-", color=color, ecolor=color, lw=1.8, ms=5.2,
            capsize=3, capthick=1.0, zorder=3,
        )
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Residual log$_{10}$ group luminosity (dex)")
        ax.grid(axis="y", color="#d8dde3", lw=0.6, alpha=0.8)
        ax.text(0.03, 0.95, f"N = {len(data)}; six equal-count bins",
                transform=ax.transAxes, va="top", color="#444444", fontsize=8.5)
    axes[0].set_ylabel("Adjusted local-expansion residual (%)")
    fig.suptitle("Binned luminosity-local-expansion trends after nuisance adjustment",
                 y=1.02, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "supplementary_figure_adjusted_binned_trends.png",
                bbox_inches="tight", dpi=300)
    plt.close(fig)


def supplementary_slope_figure() -> None:
    """Show adjusted slopes and HC3 intervals across the three analysis samples."""
    result = json.loads((RESULTS / "analysis_results.json").read_text(encoding="utf-8"))
    labels = ["Primary\n(N=858)", "Strict quality\n(N=242)", "SN Ia only\n(N=188)"]
    keys = ["primary", "strict", "sn_only"]
    beta = np.array([result[key]["adjusted_ols"]["beta"] for key in keys])
    low = np.array([result[key]["adjusted_ols"]["ci95_hc3"][0] for key in keys])
    high = np.array([result[key]["adjusted_ols"]["ci95_hc3"][1] for key in keys])

    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    x = np.arange(3)
    ax.axhline(0, color="#6f6f6f", lw=1, ls="--")
    ax.plot(x, beta, color="#315f88", lw=1.4, zorder=2)
    ax.errorbar(x, beta, yerr=[beta - low, high - beta], fmt="o", ms=7,
                color="#173f5f", ecolor="#315f88", capsize=5,
                capthick=1.2, lw=1.5, zorder=3)
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"Adjusted slope $\beta$ (dex/dex)")
    ax.set_title("Association estimate across analysis samples", fontweight="bold")
    ax.grid(axis="y", color="#d8dde3", lw=0.6, alpha=0.8)
    for index, value in enumerate(beta):
        ax.annotate(f"{value:.3f}", (index, value), xytext=(0, 10),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color="#173f5f")
    fig.tight_layout()
    fig.savefig(OUT / "supplementary_figure_sensitivity_slopes.png",
                bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    manuscript_figure()
    supplementary_binned_figure()
    supplementary_slope_figure()
    for path in sorted(OUT.glob("*.png")):
        print(path)


if __name__ == "__main__":
    main()
