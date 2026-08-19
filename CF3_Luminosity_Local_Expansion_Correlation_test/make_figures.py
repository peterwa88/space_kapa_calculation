from pathlib import Path
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cf3")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 160,
})


def binned(d, nbins=6, seed=20260819):
    d = d.copy()
    d["bin"] = pd.qcut(d["resid_logL"], nbins, labels=False, duplicates="drop")
    rng = np.random.default_rng(seed)
    rows = []
    for b, s in d.groupby("bin"):
        x = s["resid_logL"].to_numpy()
        y = s["resid_logH"].to_numpy()
        boot = []
        for _ in range(4000):
            ii = rng.integers(0, len(s), len(s))
            boot.append(y[ii].mean())
        lo, hi = np.quantile(boot, [0.025, 0.975])
        rows.append({
            "x": x.mean(),
            "y": 100 * (10 ** y.mean() - 1),
            "lo": 100 * (10 ** lo - 1),
            "hi": 100 * (10 ** hi - 1),
            "n": len(s),
        })
    return pd.DataFrame(rows)


samples = [
    ("Primary grouped-distance sample", "primary", "#1f5b99"),
    ("SN Ia-only distance sample", "sn_only", "#a65a2e"),
]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharey=True)
for ax, (title, name, color) in zip(axes, samples):
    d = pd.read_csv(ROOT / f"analysis_{name}.csv")
    q = binned(d)
    ax.axhline(0, color="#7a7a7a", lw=0.9, ls="--", zorder=1)
    ax.errorbar(q["x"], q["y"], yerr=[q["y"]-q["lo"], q["hi"]-q["y"]],
                fmt="o-", color=color, ecolor=color, lw=1.8, ms=5.2,
                capsize=3, capthick=1.0, zorder=3)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Residual log$_{10}$ group luminosity (dex)")
    ax.grid(axis="y", color="#d8dde3", lw=0.6, alpha=0.8)
    ax.text(0.03, 0.95, f"N = {len(d)}; six equal-count bins",
            transform=ax.transAxes, va="top", color="#444444", fontsize=8.5)
axes[0].set_ylabel("Adjusted local-expansion residual (%)")
fig.suptitle("Binned luminosity–local-expansion trends after nuisance adjustment",
             y=1.02, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "figure1_adjusted_binned_trends.png", bbox_inches="tight", dpi=300)
plt.close(fig)


res = json.loads((ROOT / "analysis_results.json").read_text())
labels = ["Primary\n(N=858)", "Strict quality\n(N=242)", "SN Ia only\n(N=188)"]
keys = ["primary", "strict", "sn_only"]
b = np.array([res[k]["adjusted_ols"]["beta"] for k in keys])
lo = np.array([res[k]["adjusted_ols"]["ci95_hc3"][0] for k in keys])
hi = np.array([res[k]["adjusted_ols"]["ci95_hc3"][1] for k in keys])
fig, ax = plt.subplots(figsize=(6.5, 3.4))
x = np.arange(3)
ax.axhline(0, color="#6f6f6f", lw=1, ls="--")
ax.plot(x, b, color="#315f88", lw=1.4, zorder=2)
ax.errorbar(x, b, yerr=[b-lo, hi-b], fmt="o", ms=7, color="#173f5f",
            ecolor="#315f88", capsize=5, capthick=1.2, lw=1.5, zorder=3)
ax.set_xticks(x, labels)
ax.set_ylabel("Adjusted slope β (dex/dex)")
ax.set_title("Association estimate across analysis samples", fontweight="bold")
ax.grid(axis="y", color="#d8dde3", lw=0.6, alpha=0.8)
for i, val in enumerate(b):
    ax.annotate(f"{val:.3f}", (i, val), xytext=(0, 10), textcoords="offset points",
                ha="center", fontsize=8.5, color="#173f5f")
fig.tight_layout()
fig.savefig(OUT / "figure2_sensitivity_slopes.png", bbox_inches="tight", dpi=300)
plt.close(fig)

print(OUT / "figure1_adjusted_binned_trends.png")
print(OUT / "figure2_sensitivity_slopes.png")
