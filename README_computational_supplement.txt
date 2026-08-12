Directional Modular–Weierstrass Computational Supplement
========================================================

Files
-----
1. directional_modular_weierstrass_solver.py
   High-precision solver used for the manuscript.
2. directional_modular_weierstrass_input.json
   Reproducible principal-run configuration.
3. directional_modular_weierstrass_output.json
   Complete nested machine-readable results.
4. directional_modular_weierstrass_output.csv
   Flattened results, including truncation-scan values.
5. directional_modular_weierstrass_output_summary.txt
   Human-readable benchmark summary.
6. requirements.txt
   Python dependency declaration.

Requirements
------------
Python 3.10 or newer is recommended.
mpmath version used to generate the supplied output: 1.3.0

Run
---
python directional_modular_weierstrass_solver.py \
  directional_modular_weierstrass_input.json \
  --output-json directional_modular_weierstrass_output.json \
  --output-csv directional_modular_weierstrass_output.csv

Main capabilities
-----------------
- 100-decimal-digit arbitrary-precision calculation;
- symmetric finite lattice sums for g2 and g3;
- independent Eisenstein q-series cross-check;
- real and complex near-square deformation solving;
- complex scale coordinate K and magnitude kappa = |K|;
- regional and residual branches;
- invariant-pair, discriminant, j-invariant, and closure-residual checks;
- lattice-truncation scan.

Reproducibility note
--------------------
The N=8 finite-sum value is the manuscript reproducibility benchmark. The
independent q-series result differs by approximately 0.0732%, supporting the
reported stability in order and normalization. Additional displayed digits
are computational and do not imply corresponding physical precision in the
phenomenological inputs.
