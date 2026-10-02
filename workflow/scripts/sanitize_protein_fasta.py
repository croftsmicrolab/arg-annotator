#!/usr/bin/env python3
"""Remove terminal stop symbols while preserving protein FASTA identifiers."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-empty", action="store_true", help="Allow a valid empty candidate selection")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    with args.input.open() as source, args.output.open("w") as target:
        for line in source:
            if line.startswith(">"):
                records += 1
                target.write(line)
            else:
                target.write(line.strip().rstrip("*") + "\n")
    if records == 0 and not args.allow_empty:
        raise SystemExit(f"No FASTA records were found in {args.input}")
    print(f"Sanitized {records} protein sequences")


if __name__ == "__main__":
    main()
