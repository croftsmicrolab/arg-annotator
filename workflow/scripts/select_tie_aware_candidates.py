#!/usr/bin/env python3
"""Apply the documented tie-aware RepORF selection rule and retain support-one records."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def fasta_id(header: str) -> str:
    return header[1:].split(" # ", 1)[0].split()[0]


def extract(input_path: Path, output_path: Path, selected: set[str]) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    keep = False
    with input_path.open() as source, output_path.open("w") as target:
        for line in source:
            if line.startswith(">"):
                keep = fasta_id(line) in selected
                count += int(keep)
            if keep:
                target.write(line)
    return count


def write_table(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["top_fraction", "minimum_support"], default="top_fraction")
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--fraction", required=True, type=float)
    parser.add_argument("--minimum-primary-support", required=True, type=int)
    parser.add_argument("--orf-fna", required=True, type=Path)
    parser.add_argument("--orf-faa", required=True, type=Path)
    parser.add_argument("--primary-table", required=True, type=Path)
    parser.add_argument("--primary-ids", required=True, type=Path)
    parser.add_argument("--primary-fna", required=True, type=Path)
    parser.add_argument("--primary-faa", required=True, type=Path)
    parser.add_argument("--exploratory-table", required=True, type=Path)
    parser.add_argument("--exploratory-ids", required=True, type=Path)
    parser.add_argument("--exploratory-fna", required=True, type=Path)
    parser.add_argument("--exploratory-faa", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    if not 0 < args.fraction <= 1:
        raise ValueError("--fraction must be greater than zero and no greater than one")
    with args.table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not rows:
        raise SystemExit("The complete RepORF table is empty")

    expected_order = sorted(
        rows,
        key=lambda row: (
            -int(row["reporf_support"]),
            -int(row["orf_cluster_member_count"]),
            row["reporf_id"],
        ),
    )
    if rows != expected_order:
        raise SystemExit("The RepORF table does not follow the documented deterministic ranking")

    nominal_n = max(1, math.ceil(len(rows) * args.fraction))
    boundary_support = int(rows[nominal_n - 1]["reporf_support"])
    if args.minimum_primary_support < 2:
        raise ValueError("Minimum support must be at least 2")
    if args.mode == "minimum_support":
        primary_rows = [row for row in rows if int(row["reporf_support"]) >= args.minimum_primary_support]
        strategy = "minimum_support"
    elif boundary_support > 1:
        primary_rows = rows[:nominal_n]
        strategy = "nominal_top_fraction"
    else:
        primary_rows = [
            row for row in rows if int(row["reporf_support"]) >= args.minimum_primary_support
        ]
        strategy = "support_threshold_due_to_boundary_support_one"
    exploratory_rows = [row for row in rows if int(row["reporf_support"]) == 1]

    primary_ids = {row["reporf_id"] for row in primary_rows}
    exploratory_ids = {row["reporf_id"] for row in exploratory_rows}
    write_table(args.primary_table, fields, primary_rows)
    write_table(args.exploratory_table, fields, exploratory_rows)
    args.primary_ids.parent.mkdir(parents=True, exist_ok=True)
    args.primary_ids.write_text("".join(f"{row['reporf_id']}\n" for row in primary_rows))
    args.exploratory_ids.parent.mkdir(parents=True, exist_ok=True)
    args.exploratory_ids.write_text("".join(f"{row['reporf_id']}\n" for row in exploratory_rows))

    primary_fna_count = extract(args.orf_fna, args.primary_fna, primary_ids)
    primary_faa_count = extract(args.orf_faa, args.primary_faa, primary_ids)
    exploratory_fna_count = extract(args.orf_fna, args.exploratory_fna, exploratory_ids)
    exploratory_faa_count = extract(args.orf_faa, args.exploratory_faa, exploratory_ids)
    if primary_fna_count != len(primary_rows) or primary_faa_count != len(primary_rows):
        raise SystemExit(
            "Primary sequence extraction mismatch: "
            f"table={len(primary_rows)}, fna={primary_fna_count}, faa={primary_faa_count}"
        )
    if exploratory_fna_count != len(exploratory_rows) or exploratory_faa_count != len(
        exploratory_rows
    ):
        raise SystemExit(
            "Exploratory sequence extraction mismatch: "
            f"table={len(exploratory_rows)}, fna={exploratory_fna_count}, "
            f"faa={exploratory_faa_count}"
        )

    metrics = [
        ("selection_mode", args.mode),
        ("total_reporfs", len(rows)),
        ("nominal_fraction", args.fraction),
        ("nominal_n", nominal_n),
        ("nominal_boundary_support", boundary_support),
        ("selection_strategy", strategy),
        ("minimum_primary_support", args.minimum_primary_support),
        ("primary_candidate_count", len(primary_rows)),
        ("exploratory_support1_count", len(exploratory_rows)),
    ]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)
    print(
        f"Selection strategy={strategy}; nominal_n={nominal_n}; "
        f"boundary_support={boundary_support}; primary={len(primary_rows)}; "
        f"exploratory_support1={len(exploratory_rows)}"
    )


if __name__ == "__main__":
    main()
