"""Run the complete CF3 luminosity–local-expansion reproduction pipeline."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def main() -> None:
    log_path = ROOT / "analysis_run.txt"
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(
            [sys.executable, str(ROOT / "analyze_cf3.py")],
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    subprocess.run(
        [sys.executable, str(ROOT / "make_figures.py")],
        cwd=ROOT,
        check=True,
    )
    print(f"Analysis complete. Results: {ROOT / 'analysis_results.json'}")
    print(f"Figures: {ROOT / 'figures'}")


if __name__ == "__main__":
    main()
