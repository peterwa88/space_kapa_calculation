# CF3 luminosity–local-expansion reproducibility package

This archive accompanies the manuscript:

> Peter Yongtao Wang, *An Exploratory Observational Test of a Luminosity–Local-Expansion Association in Cosmicflows-3 Galaxy Groups*.

It contains the Python analysis, the original catalog tables used by the scripts, the derived analysis tables, machine-readable statistical results, and the two manuscript figures.

## Scientific scope

The analysis defines the catalog-level local-expansion proxy

`Hloc = Vcmb / Dist`

and tests its association with the group K-s-band luminosity `logLKs`. The adjusted model controls for a cubic function of log CMB-frame velocity, the catalog luminosity correction factor, group richness, number and uncertainty of distance measurements, and a Galactic-coordinate sky dipole. HC3 standard errors, bootstrap intervals, residual Spearman statistics, and within-velocity-bin permutation tests are reported.

The primary and strict-quality group samples show a positive adjusted association. The SN Ia-only sensitivity sample is statistically consistent with zero. Accordingly, the manuscript interprets the result as preliminary catalog-level evidence, not as a causal demonstration or a precision measurement of a new physical effect.

An additional truncation sensitivity analysis repeats the primary model without the broad `40 <= Hloc <= 110` outcome cut. It retains 860 groups and tests whether the primary association was created by that cut.

## Data provenance

The two raw TSV files were downloaded without modification from the CDS VizieR catalog **J/AJ/152/50**, *Cosmicflows-3 catalog (CF3)* (Tully et al. 2016):

- `cf3_table2_all.tsv` — summary group properties (`table2`)
- `cf3_table3_all.tsv` — individual galaxy properties (`table3`)

Dataset landing page: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/152/50

Dataset DOI: https://doi.org/10.26093/cds/vizier.51520050

The TSV headers retain the exact VizieR request URLs, retrieval dates, column descriptions, citation metadata, and rights link.

## Files

- `analyze_cf3.py` — sample construction, adjusted regressions, bootstrap and permutation tests, and CSV/JSON output.
- `make_figures.py` — creates the two publication figures from the derived CSV and JSON files.
- `run_all.py` — runs the complete analysis and figure-generation pipeline.
- `requirements.txt` — Python package requirements.
- `analysis_primary.csv` — primary sample plus residualized variables (N = 858).
- `analysis_primary_no_hloc_cut.csv` — primary-eligible sample without the broad `Hloc` outcome cut (N = 860).
- `analysis_strict.csv` — strict-quality sample plus residualized variables (N = 242).
- `analysis_sn_only.csv` — SN Ia-only sensitivity sample plus residualized variables (N = 188).
- `analysis_results.json` — full machine-readable summary statistics.
- `analysis_run.txt` — console output from the verified analysis run.
- `figures/figure1_adjusted_binned_trends.png` — adjusted binned trend chart.
- `figures/figure2_sensitivity_slopes.png` — adjusted slope sensitivity chart.
- `CHECKSUMS.sha256` — SHA-256 checksums for the archived files.

## Reproduction

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

Alternatively, run the two stages separately:

```bash
python analyze_cf3.py > analysis_run.txt
python make_figures.py
```

The scripts use fixed random seed `20260819`. They overwrite the derived CSV, JSON, text, and PNG outputs in the extracted directory while leaving the raw TSV files unchanged.

## Principal verified results

| Sample | N | Adjusted beta (dex/dex) | HC3 95% CI | p-value | Approximate Hloc change per luminosity dex |
|---|---:|---:|---:|---:|---:|
| Primary | 858 | 0.0384 | [0.0068, 0.0700] | 0.0176 | +9.24% |
| Primary without broad Hloc cut | 860 | 0.0381 | [0.0061, 0.0701] | 0.0200 | +9.16% |
| Strict quality | 242 | 0.0963 | [0.0302, 0.1624] | 0.00468 | +24.83% |
| SN Ia only | 188 | 0.0106 | [-0.0159, 0.0370] | 0.435 | +2.46% |

## Citation and reuse

When reusing the raw catalog data, cite Tully et al. (2016), *The Astronomical Journal*, 152, 50, and the VizieR dataset DOI above. Consult the VizieR rights URI embedded in the TSV headers for catalog reuse conditions.

Prepared for Peter Yongtao Wang, Shandong Qizhang Intelligent Technology Co., Ltd.
