# CF3 luminosity-local-expansion-proxy correlation test

Complete reproducibility package for:

> Peter Yongtao Wang, *An Exploratory Correlation Test of Luminosity and a Local-Expansion Proxy in Cosmicflows-3 Galaxy Groups*.

Version 2.0 accompanies the refined manuscript. It contains the original catalog tables, analysis code, derived samples, machine-readable results, manuscript and supplementary figures, numerical verification, documentation, and the final Word/PDF manuscript.

## Scientific scope

The analysis defines the observer-frame catalog proxy

`H_loc = v_CMB / D`

and tests whether its nuisance-adjusted association with total group `K_s`-band luminosity has a positive sign. The proxy is derived from redshift and distance; it is not a direct covariant expansion scalar and retains peculiar-motion, calibration, selection, and catalog-construction effects. The calculation therefore tests association, not causation, a unique response function, or a value of the fixed spatial characteristic `kappa`.

## Data provenance

The two raw TSV tables in `data/raw/` were downloaded without modification from CDS VizieR catalog **J/AJ/152/50**, *Cosmicflows-3 catalog (CF3)* (Tully et al. 2016):

- `cf3_table2_all.tsv`: summary group properties (`table2`)
- `cf3_table3_all.tsv`: individual-galaxy properties (`table3`)

The embedded VizieR headers preserve the request URLs, retrieval timestamp (2026-08-18), column descriptions, citation metadata, and rights link.

- Catalog landing page: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/152/50
- Dataset DOI: https://doi.org/10.26093/cds/vizier.51520050
- Catalog article: Tully et al., *Astronomical Journal* 152, 50 (2016), https://doi.org/10.3847/0004-6256/152/2/50

## Package contents

- `code/analyze_cf3.py`: sample construction, adjusted regression, HC3 uncertainty, bootstrap, rank and permutation analyses.
- `code/make_figures.py`: final simplified manuscript figure and two supplementary diagnostic figures.
- `code/verify_outputs.py`: exact checks of the principal reported statistics.
- `code/make_checksums.py`: regenerates `CHECKSUMS.sha256`.
- `data/raw/`: original VizieR CF3 tables used by the scripts.
- `data/derived/`: primary, strict-quality, and SN Ia-only analysis tables with residualized variables.
- `results/analysis_results.json`: complete machine-readable numerical output.
- `results/analysis_run.txt`: captured output from the verified analysis.
- `figures/figure1_primary_positive_association.png`: simplified figure used in the final manuscript.
- `figures/supplementary_*.png`: additional sensitivity displays retained for auditability.
- `documentation/METHODS_AND_EQUATIONS.md`: step-by-step calculation specification.
- `documentation/DATA_DICTIONARY.md`: raw and derived variable definitions.
- `manuscript/`: final editable Word manuscript and its rendered PDF.
- `environment_verified.txt`: package versions used for the final verification.
- `CHECKSUMS.sha256`: integrity hashes for every archived file other than the checksum file itself.

## Reproduce everything

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

The fixed analysis seed is `20260819`; manuscript-figure bin intervals use seed `20260820`. `run_all.py` rebuilds the derived tables, JSON results, text log, and figures, then verifies the principal statistics.

To refresh the integrity manifest after an intentional update:

```bash
python code/make_checksums.py
```

## Principal verified results

| Sample | N | Adjusted beta (dex/dex) | HC3 95% CI | p-value | Approximate proxy change per luminosity dex |
|---|---:|---:|---:|---:|---:|
| Primary | 858 | 0.0384 | [0.0068, 0.0700] | 0.0176 | +9.24% |
| Strict quality | 242 | 0.0963 | [0.0302, 0.1624] | 0.00468 | +24.83% |
| SN Ia only | 188 | 0.0106 | [-0.0159, 0.0370] | 0.435 | +2.46% |

For the primary sample, residual Spearman `rho = 0.0850` (`p = 0.0128`) and the within-velocity-bin permutation test gives `p = 0.0144`. The SN Ia-only result is null. The defensible conclusion is preliminary sign-level consistency evidence that justifies continued testing.

## Citation and reuse

When reusing the catalog data, cite Tully et al. (2016), the VizieR dataset DOI, and follow the VizieR rights URI embedded in the TSV headers. Cite the manuscript when reusing this analysis or its derived products.

Prepared for Peter Yongtao Wang, Shandong Qizhang Intelligent Technology Co., Ltd.
