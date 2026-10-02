#!/usr/bin/env python3
"""Convert two-column MMseqs representative/member output into a size table."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    clusters: dict[str, list[str]] = defaultdict(list)
    seen: dict[str, set[str]] = defaultdict(set)
    member_owner: dict[str, str] = {}
    with args.input.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid MMseqs row at line {line_number}: {raw.rstrip()}")
            representative, member = fields[0], fields[1]
            previous = member_owner.setdefault(member, representative)
            if previous != representative:
                raise ValueError(
                    f"MMseqs member {member} occurs in multiple clusters: {previous}, {representative}"
                )
            if member not in seen[representative]:
                seen[representative].add(member)
                clusters[representative].append(member)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(clusters.items(), key=lambda item: (-len(item[1]), item[0]))
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["cluster_id", "size", "members"],
            lineterminator="\n",
        )
        writer.writeheader()
        for representative, members in ordered:
            writer.writerow(
                {"cluster_id": representative, "size": len(members), "members": ",".join(members)}
            )
    print(f"Summarized {len(ordered)} clusters containing {len(member_owner)} unique members")


if __name__ == "__main__":
    main()
