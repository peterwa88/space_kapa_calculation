"""Regenerate SHA-256 checksums for all archived files except CHECKSUMS.sha256."""

from pathlib import Path
import hashlib


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    lines = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.name != "CHECKSUMS.sha256":
            lines.append(f"{digest(path)}  {path.relative_to(ROOT).as_posix()}")
    (ROOT / "CHECKSUMS.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(lines)} checksums.")


if __name__ == "__main__":
    main()
