# Data dictionary

## Raw CF3 fields used

| Field | Source | Meaning in this analysis |
|---|---|---|
| `Nest` | tables 2 and 3 | CF3/2MASS group identifier |
| `<Dist>` | table 2 | Group luminosity distance in Mpc; renamed `Dist` |
| `<Vcmb>` | tables 2 and 3 | CMB-frame group velocity in km/s; renamed `Vcmb` |
| `logLKs` | tables 2 and 3 | Base-10 logarithm of summed group `K_s` luminosity in solar units |
| `CF` | tables 2 and 3 | Catalog luminosity selection-function correction factor |
| `Npv` | tables 2 and 3 | Number of group galaxies with positions and velocities |
| `o_<DM>` | tables 2 and 3 | Number of contributing distance measurements |
| `e_<DM>` | tables 2 and 3 | Group distance-modulus uncertainty in magnitudes |
| `GLON`, `GLAT` | table 2 | Galactic longitude and latitude in degrees |
| `GGLON`, `GGLAT` | table 3 | Group Galactic longitude and latitude in degrees |
| `DM-N` | table 3 | SN Ia distance modulus used for the SN-only sensitivity sample |

The raw files retain all downloaded VizieR columns, even though the computation uses only the fields above.

## Derived variables

| Variable | Definition |
|---|---|
| `Hloc` | `Vcmb / Dist` |
| `logH` | `log10(Hloc)` |
| `logV` | `log10(Vcmb)` |
| `logCF` | `log10(CF)` |
| `logNpv` | `log10(Npv)` |
| `logNdist` | `log10(number of distance measurements)` |
| `eDM` | Distance-modulus uncertainty |
| `sky_x` | `cos(GLAT) cos(GLON)` |
| `sky_y` | `cos(GLAT) sin(GLON)` |
| `sky_z` | `sin(GLAT)` |
| `resid_logL` | Residual of `logLKs` after the full nuisance projection |
| `resid_logH` | Residual of `logH` after the same nuisance projection |

Angles are converted to radians before the sky-vector calculation.

## Derived tables

- `analysis_primary.csv`: 858 primary-sample groups.
- `analysis_strict.csv`: 242 nested strict-quality groups.
- `analysis_sn_only.csv`: 188 SN Ia-only group aggregates.

Each derived table contains the retained raw columns, constructed variables, and the two nuisance-adjusted residual columns required to reproduce the figures.
