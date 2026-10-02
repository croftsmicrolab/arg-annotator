#!/usr/bin/env python3
"""Compute RepORF support from unique parent insert clusters."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def fasta_headers(path: Path):
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].rstrip("\n")


def parse_orf_headers(path: Path):
    mapping: dict[str, str] = {}
    metadata: dict[str, dict[str, str]] = {}
    for header in fasta_headers(path):
        parts = [part.strip() for part in header.split(" # ")]
        orf_id = parts[0].split()[0]
        match = re.match(r"^(.+)_([0-9]+)$", orf_id)
        if not match:
            raise ValueError(f"Cannot infer the parent insert from Prodigal ORF ID: {orf_id}")
        if orf_id in mapping:
            raise ValueError(f"Duplicate ORF ID in FASTA: {orf_id}")
        parent = match.group(1)
        attributes = parts[4] if len(parts) > 4 else ""
        partial_match = re.search(r"(?:^|;)partial=([^;]+)", attributes)
        mapping[orf_id] = parent
        metadata[orf_id] = {
            "start": parts[1] if len(parts) > 1 else "",
            "end": parts[2] if len(parts) > 2 else "",
            "strand": parts[3] if len(parts) > 3 else "",
            "partial_flag": partial_match.group(1) if partial_match else "",
        }
    return mapping, metadata


def parse_insert_summary(path: Path):
    member_to_cluster: dict[str, str] = {}
    sizes: dict[str, int] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cluster = row["cluster_id"]
            size = int(row["size"])
            members = [member for member in row["members"].split(",") if member]
            if size != len(members):
                raise ValueError(
                    f"Insert cluster size mismatch for {cluster}: table={size}, members={len(members)}"
                )
            sizes[cluster] = size
            for member in members:
                previous = member_to_cluster.setdefault(member, cluster)
                if previous != cluster:
                    raise ValueError(f"Insert {member} occurs in multiple clusters")
    return member_to_cluster, sizes


def parse_orf_clusters(path: Path):
    clusters: dict[str, list[str]] = defaultdict(list)
    seen_by_cluster: dict[str, set[str]] = defaultdict(set)
    member_owner: dict[str, str] = {}
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid ORF-cluster row at line {line_number}")
            representative, member = fields[0], fields[1]
            previous = member_owner.setdefault(member, representative)
            if previous != representative:
                raise ValueError(f"ORF {member} occurs in multiple clusters")
            if member not in seen_by_cluster[representative]:
                seen_by_cluster[representative].add(member)
                clusters[representative].append(member)
    return clusters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orfs", required=True, type=Path)
    parser.add_argument("--orf-clusters", required=True, type=Path)
    parser.add_argument("--insert-clusters", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mapping-output", required=True, type=Path)
    parser.add_argument("--missing-output", required=True, type=Path)
    args = parser.parse_args()

    orf_to_parent, orf_metadata = parse_orf_headers(args.orfs)
    member_to_insert_cluster, insert_sizes = parse_insert_summary(args.insert_clusters)
    orf_clusters = parse_orf_clusters(args.orf_clusters)

    for path in (args.output, args.mapping_output, args.missing_output):
        path.parent.mkdir(parents=True, exist_ok=True)

    mapping_fields = [
        "orf_id",
        "parent_insert_rep",
        "insert_cluster_id",
        "insert_cluster_size",
        "start",
        "end",
        "strand",
        "partial_flag",
    ]
    with args.mapping_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=mapping_fields, lineterminator="\n")
        writer.writeheader()
        for orf_id in sorted(orf_to_parent):
            parent = orf_to_parent[orf_id]
            cluster = member_to_insert_cluster.get(parent, "")
            writer.writerow(
                {
                    "orf_id": orf_id,
                    "parent_insert_rep": parent,
                    "insert_cluster_id": cluster,
                    "insert_cluster_size": insert_sizes.get(cluster, ""),
                    **orf_metadata[orf_id],
                }
            )

    output_rows = []
    missing_rows = []
    for representative, members in orf_clusters.items():
        parents: set[str] = set()
        insert_clusters: set[str] = set()
        for orf_id in members:
            parent = orf_to_parent.get(orf_id)
            if parent is None:
                missing_rows.append(("orf_not_found_in_fasta", representative, orf_id, ""))
                continue
            parents.add(parent)
            cluster = member_to_insert_cluster.get(parent)
            if cluster is None:
                missing_rows.append(
                    ("parent_not_found_in_insert_clusters", representative, orf_id, parent)
                )
                continue
            insert_clusters.add(cluster)
        support = sum(insert_sizes[cluster] for cluster in insert_clusters)
        output_rows.append(
            {
                "reporf_id": representative,
                "orf_cluster_member_count": len(members),
                "unique_parent_insert_count": len(parents),
                "unique_insert_cluster_count": len(insert_clusters),
                "reporf_support": support,
                "parent_insert_reps": ",".join(sorted(parents)),
                "insert_cluster_ids": ",".join(sorted(insert_clusters)),
                "member_orfs": ",".join(sorted(members)),
            }
        )

    output_rows.sort(
        key=lambda row: (
            -int(row["reporf_support"]),
            -int(row["orf_cluster_member_count"]),
            str(row["reporf_id"]),
        )
    )
    support_fields = [
        "reporf_rank",
        "reporf_id",
        "orf_cluster_member_count",
        "unique_parent_insert_count",
        "unique_insert_cluster_count",
        "reporf_support",
        "parent_insert_reps",
        "insert_cluster_ids",
        "member_orfs",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=support_fields, lineterminator="\n")
        writer.writeheader()
        for rank, row in enumerate(output_rows, 1):
            writer.writerow({"reporf_rank": rank, **row})

    with args.missing_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["problem", "reporf_id", "orf_id", "parent_insert_rep"])
        writer.writerows(missing_rows)

    print(f"Ranked {len(output_rows)} RepORFs; missing_mappings={len(missing_rows)}")
    if missing_rows:
        raise SystemExit("Incomplete ORF-to-insert mapping; inspect missing_mappings.tsv")


if __name__ == "__main__":
    main()
