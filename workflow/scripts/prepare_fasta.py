#!/usr/bin/env python3
"""Normalize FASTA identifiers and write record-level nucleotide QC."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def records(path: Path):
    header = None
    sequence = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence encountered before the first FASTA header in {path}")
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def clean_id(original: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", original).strip("_")
    return value or "sequence"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--allowed-symbols", required=True)
    parser.add_argument("--line-width", required=True, type=int)
    parser.add_argument("--fail-on-invalid-bases", action="store_true")
    parser.add_argument("--fail-on-empty-sequences", action="store_true")
    args = parser.parse_args()

    if args.line_width < 1:
        raise ValueError("--line-width must be at least 1")
    allowed = set(args.allowed_symbols.upper())
    for path in (args.output, args.mapping, args.qc):
        path.parent.mkdir(parents=True, exist_ok=True)

    seen: dict[str, int] = {}
    used_ids = set()
    count = 0
    invalid_total = 0
    empty_total = 0
    with (
        args.output.open("w") as fasta_out,
        args.mapping.open("w", newline="") as map_out,
        args.qc.open("w", newline="") as qc_out,
    ):
        map_writer = csv.DictWriter(
            map_out,
            delimiter="\t",
            fieldnames=["original_sequence_id", "cleaned_sequence_id", "original_description"],
            lineterminator="\n",
        )
        qc_writer = csv.DictWriter(
            qc_out,
            delimiter="\t",
            fieldnames=[
                "sequence_id",
                "sequence_length",
                "gc_content",
                "n_count",
                "n_fraction",
                "invalid_base_count",
            ],
            lineterminator="\n",
        )
        map_writer.writeheader()
        qc_writer.writeheader()

        for description, raw_sequence in records(args.input):
            original_id = description.split()[0]
            base = clean_id(original_id)
            seen[base] = seen.get(base, 0) + 1
            cleaned = base if seen[base] == 1 else f"{base}__dup{seen[base]}"
            while cleaned in used_ids:
                seen[base] += 1
                cleaned = f"{base}__dup{seen[base]}"
            used_ids.add(cleaned)
            sequence = raw_sequence.upper()
            invalid = sum(symbol not in allowed for symbol in sequence)
            gc_count = sequence.count("G") + sequence.count("C")
            n_count = sequence.count("N")
            length = len(sequence)
            invalid_total += invalid
            empty_total += int(length == 0)

            fasta_out.write(f">{cleaned}\n")
            for start in range(0, length, args.line_width):
                fasta_out.write(sequence[start : start + args.line_width] + "\n")
            map_writer.writerow(
                {
                    "original_sequence_id": original_id,
                    "cleaned_sequence_id": cleaned,
                    "original_description": description,
                }
            )
            qc_writer.writerow(
                {
                    "sequence_id": cleaned,
                    "sequence_length": length,
                    "gc_content": f"{(100 * gc_count / length) if length else 0:.6f}",
                    "n_count": n_count,
                    "n_fraction": f"{(n_count / length) if length else 0:.8f}",
                    "invalid_base_count": invalid,
                }
            )
            count += 1

    if count == 0:
        raise SystemExit(f"No FASTA records were found in {args.input}")
    if args.fail_on_empty_sequences and empty_total:
        raise SystemExit(f"QC failed: {empty_total} empty sequences were found")
    if args.fail_on_invalid_bases and invalid_total:
        raise SystemExit(f"QC failed: {invalid_total} invalid nucleotide symbols were found")
    print(
        f"Prepared {count} sequences; empty_sequences={empty_total}; "
        f"invalid_base_count={invalid_total}"
    )


if __name__ == "__main__":
    main()
