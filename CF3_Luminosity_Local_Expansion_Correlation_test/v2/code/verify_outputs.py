"""Verify the sample sizes and principal manuscript statistics."""

from pathlib import Path
import json
import math


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "analysis_results.json"

EXPECTED = {
    "primary": {"n": 858, "beta": 0.038366042231306664, "p": 0.01755080082602992},
    "strict": {"n": 242, "beta": 0.09633212533444965, "p": 0.004675581877553083},
    "sn_only": {"n": 188, "beta": 0.010563586630220556, "p": 0.4345601603869794},
}


def main() -> None:
    result = json.loads(RESULTS.read_text(encoding="utf-8"))
    for sample, expected in EXPECTED.items():
        observed = result[sample]
        assert observed["n"] == expected["n"], (sample, "n", observed["n"])
        assert math.isclose(observed["adjusted_ols"]["beta"], expected["beta"],
                            rel_tol=0.0, abs_tol=1e-12), (sample, "beta")
        assert math.isclose(observed["adjusted_ols"]["p"], expected["p"],
                            rel_tol=0.0, abs_tol=1e-12), (sample, "p")
    assert math.isclose(result["primary"]["residual_spearman"]["rho"],
                        0.08497102472621482, rel_tol=0.0, abs_tol=1e-12)
    print("Verification passed: all principal statistics match the manuscript.")


if __name__ == "__main__":
    main()
