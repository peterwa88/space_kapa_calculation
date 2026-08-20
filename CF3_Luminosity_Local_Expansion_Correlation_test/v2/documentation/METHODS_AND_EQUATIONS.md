# Calculation methods and equations

This file states the complete calculation implemented by `code/analyze_cf3.py`. It is an audit specification, not an additional theoretical claim.

## 1. Catalog proxy

For group `i`,

\[
H_{\mathrm{loc},i}=\frac{v_{\mathrm{CMB},i}}{D_i},
\qquad
y_i=\log_{10}H_{\mathrm{loc},i},
\qquad
x_i=\log_{10}L_{K_s,i}.
\]

`v_CMB` is the catalog CMB-frame recession velocity and `D` is the catalog group distance. `H_loc` is an observer-frame redshift-distance proxy, not a direct covariant expansion scalar.

## 2. Primary sample construction

Rows must contain finite `Nest`, `D`, `v_CMB`, `logLKs`, `CF`, `Npv`, distance-modulus uncertainty, and number of distance measurements. The prespecified filters are:

\[
3000\leq v_{\mathrm{CMB}}\leq15000\ \mathrm{km\,s^{-1}},
\]

\[
N_{\mathrm{dist}}\geq2,
\qquad
\sigma_\mu\leq0.30\ \mathrm{mag},
\qquad
1\leq C_F\leq10,
\]

\[
40\leq H_{\mathrm{loc}}\leq110\ \mathrm{km\,s^{-1}\,Mpc^{-1}}.
\]

These filters retain 858 groups. The final bound removes catastrophic velocity-distance ratios; it is not fitted to maximize the luminosity coefficient.

The strict sample additionally requires

\[
N_{\mathrm{dist}}\geq3,
\qquad
\sigma_\mu\leq0.20\ \mathrm{mag},
\qquad
C_F\leq5,
\]

retaining 242 groups.

The SN Ia-only sample is built from table 3 by taking the median repeated SN Ia distance modulus within each `Nest`, converting it through

\[
D=10^{(\mu-25)/5}\ \mathrm{Mpc},
\]

and applying the velocity, correction-factor, and proxy bounds above. It retains 188 groups.

## 3. Nuisance-adjusted model

The primary regression is

\[
y_i=\alpha+\beta x_i+f_3(\log_{10}v_{\mathrm{CMB},i})
+\gamma_1\log_{10}C_{F,i}
+\gamma_2\log_{10}N_{\mathrm{pv},i}
+\gamma_3\log_{10}N_{\mathrm{dist},i}
+\gamma_4\sigma_{\mu,i}
+\boldsymbol{\delta}^{\mathsf T}\mathbf{s}_i+\varepsilon_i,
\]

where `f_3` is a centered cubic polynomial and the Galactic sky-direction vector is

\[
\mathbf{s}_i=
(\cos b_i\cos l_i,\ \cos b_i\sin l_i,\ \sin b_i).
\]

The coefficient `beta` is the conditional slope in dex of `H_loc` per dex of cataloged luminosity. It is not interpreted as a causal fraction.

## 4. HC3 covariance and interval

Ordinary least squares gives

\[
\widehat{\boldsymbol\theta}=(X^{\mathsf T}X)^{-1}X^{\mathsf T}y.
\]

With residual `e_i` and leverage `h_ii`, the HC3 sandwich covariance used by the code is

\[
\widehat{\mathrm{Var}}_{\mathrm{HC3}}(\widehat{\boldsymbol\theta})
=(X^{\mathsf T}X)^{-1}
X^{\mathsf T}\operatorname{diag}\!\left[\left(\frac{e_i}{1-h_{ii}}\right)^2\right]X
(X^{\mathsf T}X)^{-1}.
\]

The reported two-sided `p` value uses the Student `t` distribution with `n-p` degrees of freedom. The displayed 95% interval is `beta +/- 1.96 SE_HC3`, matching the manuscript calculation.

The back-transformed percent change for one luminosity decade is

\[
\Delta H_{\mathrm{loc}}=100(10^{\beta}-1)\%.
\]

## 5. Rank, residual, bootstrap, and permutation checks

1. Raw Pearson and Spearman statistics are calculated from `x` and `y`.
2. Partial Spearman correlation rank-residualizes `x` and `y` against ranked `logV` and `logCF`.
3. Full-model residuals are obtained by regressing both `x` and `y` on the same nuisance matrix. Their Spearman correlation is the residual-rank check.
4. The adjusted slope is bootstrapped with 3000 row-resampling replicates and fixed seed `20260819`.
5. The permutation test shuffles residual luminosity within 20 equal-count `logV` bins, uses 5000 permutations and the same seed, and reports

\[
p_{\mathrm{perm}}=\frac{1+\#\{|\rho^{(b)}|\geq|\rho_{\mathrm{obs}}|\}}{B+1}.
\]

## 6. Figure calculation

For the final manuscript figure, primary-sample nuisance-adjusted luminosity residuals are split into six equal-count bins. Each point is the mean residual pair, with the vertical value back-transformed to percent. Error bars are percentile 95% intervals from 4000 bootstrap resamples per bin. The straight line is the Frisch-Waugh-Lovell residual regression and is a visualization of the prespecified adjusted model, not a separately selected functional law.

## 7. Interpretation boundary

The calculation supports only a weak positive conditional association in the primary CF3 sample. It does not isolate metric expansion from peculiar motion, prove luminosity-driven causation, determine a unique macroscopic response function, estimate `kappa`, or test a pulsar-timing scalar breathing mode. Those are separate future tests.
