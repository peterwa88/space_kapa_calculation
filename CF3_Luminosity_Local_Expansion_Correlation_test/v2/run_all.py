"""Run the complete CF3 analysis, figure generation, and numerical checks."""

from pathlib import Path
import os
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
ENV = os.environ.copy()
ENV.setdefault("MPLCONFIGDIR", "/tmp/mpl-cf3-v2")


def main() -> None:
    log_path = ROOT / "results" / "analysis_run.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(
            [sys.executable, str(ROOT / "code" / "analyze_cf3.py")],
            cwd=ROOT,
            env=ENV,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    subprocess.run(
        [sys.executable, str(ROOT / "code" / "make_figures.py")],
        cwd=ROOT,
        env=ENV,
        check=True,
    )
    subprocess.run(
        [sys.executable, str(ROOT / "code" / "verify_outputs.py")],
        cwd=ROOT,
        env=ENV,
        check=True,
    )
    print(f"Analysis complete: {ROOT / 'results' / 'analysis_results.json'}")
    print(f"Figures complete: {ROOT / 'figures'}")


if __name__ == "__main__":
    main()
