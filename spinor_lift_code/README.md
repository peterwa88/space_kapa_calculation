# Spinor-Lift Kappa Estimate — Reproducibility Package

This package reproduces the numerical estimates in the accompanying English Word note.

## Hypothesis evaluated

The proposed scalar spinor-lift ansatz is

`P'^mu = kappa^2 Lambda^mu_nu P^nu`.

Because a Lorentz transformation preserves the Minkowski norm, a timelike state pair that is *physically* related by the additional scalar lift would satisfy

`kappa = sqrt(M/m)`.

This is a conditional estimator, not a standard result of spinor-helicity theory. The Planck and reduced-Planck masses are used only as benchmark candidate partner scales; no known experiment has established such a partner state.

## Run

```bash
python3 spinor_lift_kappa_estimate.py
```

The script uses only the Python standard library and writes:

- `spinor_lift_results.json`
- `spinor_lift_results.csv`

## Core numerical inputs

- Planck mass energy equivalent: `1.220890e19 GeV` (2022 CODATA/NIST central value).
- Atmospheric mass-squared splitting benchmark: `2.513e-3 eV^2` (NuFIT 6.0, normal-ordering benchmark).
- KATRIN direct effective-neutrino-mass upper limit: `m_beta < 0.45 eV` at 90% C.L.
- Companion modular-closure benchmark: `kappa = 2.81005508341944e14`.

## Key outputs

For `m = 0.05 eV`:

- ordinary Planck benchmark: `kappa ≈ 4.9414e14`
- reduced-Planck benchmark: `kappa ≈ 2.20695e14`

Using `kappa = 2.81005508341944e14` inversely:

- ordinary Planck benchmark predicts `m ≈ 0.15461 eV`
- reduced-Planck benchmark predicts `m ≈ 0.03084 eV`

The numerical proximity of the reduced-Planck benchmark to the modular value is not evidence by itself; the physical identification of a spinor partner scale must be derived independently.
