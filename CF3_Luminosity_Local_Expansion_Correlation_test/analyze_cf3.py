from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent


def read_vizier(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    first = df.columns[0]
    df[first] = pd.to_numeric(df[first], errors="coerce")
    return df.loc[df[first].notna()].copy()


def num(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def partial_spearman(y, x, controls):
    ry = stats.rankdata(y)
    rx = stats.rankdata(x)
    z = np.column_stack([stats.rankdata(controls[:, j]) for j in range(controls.shape[1])])
    z = np.column_stack([np.ones(len(z)), z])
    ey = ry - z @ np.linalg.lstsq(z, ry, rcond=None)[0]
    ex = rx - z @ np.linalg.lstsq(z, rx, rcond=None)[0]
    return stats.pearsonr(ex, ey)


def ols_hc3(y, X):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    inv = np.linalg.inv(X.T @ X)
    resid = y - X @ beta
    lev = np.einsum("ij,jk,ik->i", X, inv, X)
    u = resid / np.clip(1 - lev, 1e-8, None)
    meat = X.T @ ((u * u)[:, None] * X)
    cov = inv @ meat @ inv
    se = np.sqrt(np.diag(cov))
    df = len(y) - X.shape[1]
    tval = beta / se
    pval = 2 * stats.t.sf(np.abs(tval), df)
    return beta, se, pval, resid


def within_bin_permutation(y, x, bin_var, nperm=10000, seed=20260819):
    rng = np.random.default_rng(seed)
    bins = pd.qcut(bin_var, 20, labels=False, duplicates="drop").to_numpy()
    obs = stats.spearmanr(x, y).statistic
    vals = np.empty(nperm)
    for k in range(nperm):
        xp = x.copy()
        for b in np.unique(bins):
            idx = np.where(bins == b)[0]
            xp[idx] = rng.permutation(xp[idx])
        vals[k] = stats.spearmanr(xp, y).statistic
    p = (1 + np.sum(np.abs(vals) >= abs(obs))) / (nperm + 1)
    return obs, p, vals


def summarize(df, label):
    y = df["logH"].to_numpy()
    x = df["logLKs"].to_numpy()
    pear = stats.pearsonr(x, y)
    spear = stats.spearmanr(x, y)
    ps = partial_spearman(y, x, df[["logV", "logCF"]].to_numpy())

    # Flexible control for distance/redshift selection: cubic B-spline of logV,
    # plus catalog selection correction, richness, and distance-modulus error.
    lv = df["logV"].to_numpy()
    lv0 = lv - lv.mean()
    nuisance_cols = [np.ones(len(df)), lv0, lv0**2, lv0**3,
                     df["logCF"].to_numpy(), df["logNpv"].to_numpy()]
    for c in ["logNdist", "eDM", "sky_x", "sky_y", "sky_z"]:
        if c in df and df[c].std() > 1e-10:
            nuisance_cols.append(df[c].to_numpy())
    nuisance = np.column_stack(nuisance_cols)
    X = np.column_stack([np.ones(len(df)), x, nuisance[:, 1:]])
    beta, ses, pvals, _ = ols_hc3(y, X)
    b = beta[1]
    se = ses[1]
    p = pvals[1]
    ci = [b - 1.96 * se, b + 1.96 * se]

    # Residualized correlation and within-redshift-bin permutation.
    yres = y - nuisance @ np.linalg.lstsq(nuisance, y, rcond=None)[0]
    xres = x - nuisance @ np.linalg.lstsq(nuisance, x, rcond=None)[0]
    rres = stats.spearmanr(xres, yres)
    perm_obs, perm_p, _ = within_bin_permutation(yres.copy(), xres.copy(), df["logV"], nperm=5000)

    # Bootstrap the adjusted coefficient.
    rng = np.random.default_rng(20260819)
    boots = []
    for _ in range(3000):
        ii = rng.integers(0, len(df), len(df))
        try:
            boots.append(np.linalg.lstsq(X[ii], y[ii], rcond=None)[0][1])
        except Exception:
            continue
    bci = np.quantile(boots, [0.025, 0.975]).tolist()

    # Distance/richness strata for sign stability.
    strata = []
    for name, sub in df.groupby(pd.qcut(df["logV"], 4, labels=["Q1", "Q2", "Q3", "Q4"])):
        rr = stats.spearmanr(sub["logLKs"], sub["logH"])
        strata.append({"stratum": str(name), "n": len(sub), "rho": rr.statistic, "p": rr.pvalue,
                       "v_min": sub["Vcmb"].min(), "v_max": sub["Vcmb"].max()})

    return {
        "label": label,
        "n": len(df),
        "H_median": float(df["Hloc"].median()),
        "H_iqr": [float(df["Hloc"].quantile(.25)), float(df["Hloc"].quantile(.75))],
        "logL_range": [float(df["logLKs"].min()), float(df["logLKs"].max())],
        "V_range": [float(df["Vcmb"].min()), float(df["Vcmb"].max())],
        "raw_pearson": {"r": pear.statistic, "p": pear.pvalue},
        "raw_spearman": {"rho": spear.statistic, "p": spear.pvalue},
        "partial_spearman_logV_logCF": {"rho": ps.statistic, "p": ps.pvalue},
        "adjusted_ols": {"beta": b, "se_hc3": se, "p": p, "ci95_hc3": ci, "bootstrap_ci95": bci,
                         "interpret_pct_H_per_dex": 100*(10**b-1)},
        "residual_spearman": {"rho": rres.statistic, "p": rres.pvalue},
        "permutation": {"rho": perm_obs, "p": perm_p, "nperm": 5000},
        "quartiles": strata,
    }, xres, yres


g = read_vizier(ROOT / "cf3_table2_all.tsv")
num(g, ["Nest", "o_<DM>", "<DM>", "e_<DM>", "<Dist>", "Npv", "logLKs", "CF", "<Vcmb>", "GLON", "GLAT"])
g = g.rename(columns={"<Dist>": "Dist", "<Vcmb>": "Vcmb", "e_<DM>": "eDM", "o_<DM>": "Ndist"})
g["Hloc"] = g["Vcmb"] / g["Dist"]
g["logH"] = np.log10(g["Hloc"])
g["logV"] = np.log10(g["Vcmb"])
g["logCF"] = np.log10(g["CF"])
g["logNpv"] = np.log10(g["Npv"])
g["logNdist"] = np.log10(g["Ndist"])
lon = np.deg2rad(g["GLON"]); lat = np.deg2rad(g["GLAT"])
g["sky_x"] = np.cos(lat) * np.cos(lon)
g["sky_y"] = np.cos(lat) * np.sin(lon)
g["sky_z"] = np.sin(lat)

# Preregistered-style quality cuts: Hubble-flow velocities suppress peculiar-motion
# dominance; repeated distances and modulus-error cuts suppress noisy groups; CF cut
# limits the largest catalog lost-light corrections.  The broad H cut in the primary
# sample removes catastrophes; primary_no_hloc_cut retains the same eligible groups
# without that outcome-variable cut as a prespecified sensitivity analysis.
primary_no_hloc_cut = g.loc[
    g[["Nest", "Dist", "Vcmb", "logLKs", "CF", "Npv", "eDM", "Ndist"]].notna().all(axis=1)
    & g["Vcmb"].between(3000, 15000)
    & (g["Ndist"] >= 2)
    & (g["eDM"] <= 0.30)
    & g["CF"].between(1, 10)
].copy()
primary = primary_no_hloc_cut.loc[primary_no_hloc_cut["Hloc"].between(40, 110)].copy()

# A stricter high-quality sample.
strict = primary.loc[(primary["Ndist"] >= 3) & (primary["eDM"] <= 0.20) & (primary["CF"] <= 5)].copy()

# SN-Ia-only distance sensitivity sample from table 3; collapse repeated rows by nest.
t = read_vizier(ROOT / "cf3_table3_all.tsv")
num(t, ["Nest", "DM-N", "o_SN", "logLKs", "CF", "Npv", "<Vcmb>", "o_<DM>", "e_<DM>", "GGLON", "GGLAT"])
t = t.loc[t[["Nest", "DM-N", "logLKs", "CF", "Npv", "<Vcmb>"]].notna().all(axis=1)].copy()
sn = t.groupby("Nest", as_index=False).agg({
    "DM-N": "median", "logLKs": "first", "CF": "first", "Npv": "first", "<Vcmb>": "first",
    "o_<DM>": "first", "e_<DM>": "first", "GGLON": "first", "GGLAT": "first"
})
sn = sn.rename(columns={"<Vcmb>": "Vcmb"})
sn["Dist"] = 10 ** ((sn["DM-N"] - 25) / 5)
sn["Hloc"] = sn["Vcmb"] / sn["Dist"]
sn["logH"] = np.log10(sn["Hloc"])
sn["logV"] = np.log10(sn["Vcmb"])
sn["logCF"] = np.log10(sn["CF"])
sn["logNpv"] = np.log10(sn["Npv"])
sn["logNdist"] = np.log10(sn["o_<DM>"])
sn["eDM"] = sn["e_<DM>"]
lon = np.deg2rad(sn["GGLON"]); lat = np.deg2rad(sn["GGLAT"])
sn["sky_x"] = np.cos(lat) * np.cos(lon)
sn["sky_y"] = np.cos(lat) * np.sin(lon)
sn["sky_z"] = np.sin(lat)
sn = sn.loc[sn["Vcmb"].between(3000, 15000) & sn["CF"].between(1, 10) & sn["Hloc"].between(40, 110)].copy()

results = {}
samples = [
    ("primary", primary),
    ("strict", strict),
    ("sn_only", sn),
    ("primary_no_hloc_cut", primary_no_hloc_cut),
]
for label, df in samples:
    res, xr, yr = summarize(df.reset_index(drop=True), label)
    results[label] = res
    out = df.reset_index(drop=True).copy()
    out["resid_logL"] = xr
    out["resid_logH"] = yr
    out.to_csv(ROOT / f"analysis_{label}.csv", index=False)

(ROOT / "analysis_results.json").write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(json.dumps(results, indent=2, default=float))
