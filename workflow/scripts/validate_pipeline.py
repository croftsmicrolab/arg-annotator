#!/usr/bin/env python3
"""Validate sequence, clustering, support, selection, and annotation invariants."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def fasta_ids(path: Path) -> list[str]:
    return [
        line[1:].split(" # ", 1)[0].split()[0]
        for line in path.open()
        if line.startswith(">")
    ]


def table_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cluster_pairs(path: Path):
    pairs = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid cluster membership row at line {line_number}: {path}")
            pairs.append((fields[0], fields[1]))
    return pairs


def selection_metrics(path: Path) -> dict[str, str]:
    return {row["metric"]: row["value"] for row in table_rows(path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--normalized-fasta", required=True, type=Path)
    parser.add_argument("--id-mapping", required=True, type=Path)
    parser.add_argument("--sequence-qc", required=True, type=Path)
    parser.add_argument("--insert-members", required=True, type=Path)
    parser.add_argument("--insert-summary", required=True, type=Path)
    parser.add_argument("--insert-representatives", required=True, type=Path)
    parser.add_argument("--orf-fna", required=True, type=Path)
    parser.add_argument("--orf-faa", required=True, type=Path)
    parser.add_argument("--orf-members", required=True, type=Path)
    parser.add_argument("--orf-summary", required=True, type=Path)
    parser.add_argument("--representative-orfs", required=True, type=Path)
    parser.add_argument("--orf-mapping", required=True, type=Path)
    parser.add_argument("--missing-mappings", required=True, type=Path)
    parser.add_argument("--support-table", required=True, type=Path)
    parser.add_argument("--selection-summary", required=True, type=Path)
    parser.add_argument("--primary-table", required=True, type=Path)
    parser.add_argument("--primary-fna", required=True, type=Path)
    parser.add_argument("--primary-faa", required=True, type=Path)
    parser.add_argument("--exploratory-table", required=True, type=Path)
    parser.add_argument("--exploratory-fna", required=True, type=Path)
    parser.add_argument("--exploratory-faa", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--expected-input", required=True, type=int)
    parser.add_argument("--expected-insert-clusters", required=True, type=int)
    parser.add_argument("--expected-orfs", required=True, type=int)
    parser.add_argument("--expected-reporfs", required=True, type=int)
    parser.add_argument("--expected-primary", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    checks: list[tuple[str, str, str]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append((name, "PASS" if condition else "FAIL", detail))

    raw_ids = fasta_ids(args.input_fasta)
    normalized_ids = fasta_ids(args.normalized_fasta)
    mapping_rows = table_rows(args.id_mapping)
    qc_rows = table_rows(args.sequence_qc)
    mapping_ids = [row["cleaned_sequence_id"] for row in mapping_rows]
    qc_ids = [row["sequence_id"] for row in qc_rows]
    check(
        "input_sequence_count",
        len(raw_ids) == len(normalized_ids) and (args.expected_input == 0 or len(raw_ids) == args.expected_input),
        f"raw={len(raw_ids)}; normalized={len(normalized_ids)}; expected={args.expected_input}",
    )
    check(
        "input_id_normalization_complete",
        normalized_ids == mapping_ids == qc_ids and len(normalized_ids) == len(set(normalized_ids)),
        f"normalized={len(normalized_ids)}; mapping={len(mapping_ids)}; qc={len(qc_ids)}",
    )
    invalid_count = sum(int(row["invalid_base_count"]) for row in qc_rows)
    empty_count = sum(int(row["sequence_length"]) == 0 for row in qc_rows)
    check(
        "input_fasta_qc",
        invalid_count == 0 and empty_count == 0,
        f"invalid_base_count={invalid_count}; empty_sequences={empty_count}",
    )

    insert_pairs = cluster_pairs(args.insert_members)
    insert_member_counts = Counter(member for _, member in insert_pairs)
    insert_clusters_from_pairs: dict[str, list[str]] = defaultdict(list)
    for representative, member in insert_pairs:
        insert_clusters_from_pairs[representative].append(member)
    check(
        "insert_cluster_membership_complete",
        set(insert_member_counts) == set(normalized_ids)
        and all(count == 1 for count in insert_member_counts.values()),
        f"input={len(normalized_ids)}; unique_members={len(insert_member_counts)}; rows={len(insert_pairs)}",
    )

    insert_summary_rows = table_rows(args.insert_summary)
    insert_summary = {row["cluster_id"]: row for row in insert_summary_rows}
    insert_summary_members = {
        cluster: [member for member in row["members"].split(",") if member]
        for cluster, row in insert_summary.items()
    }
    insert_summary_consistent = set(insert_summary) == set(insert_clusters_from_pairs)
    if insert_summary_consistent:
        for cluster in insert_summary:
            insert_summary_consistent &= int(insert_summary[cluster]["size"]) == len(
                insert_summary_members[cluster]
            )
            insert_summary_consistent &= set(insert_summary_members[cluster]) == set(
                insert_clusters_from_pairs[cluster]
            )
    check(
        "insert_cluster_summary_consistent",
        insert_summary_consistent and (args.expected_insert_clusters == 0 or len(insert_summary) == args.expected_insert_clusters),
        f"observed_clusters={len(insert_summary)}; expected_clusters={args.expected_insert_clusters}",
    )
    insert_rep_ids = fasta_ids(args.insert_representatives)
    check(
        "insert_representatives_complete",
        set(insert_rep_ids) == set(insert_summary) and len(insert_rep_ids) == len(insert_summary),
        f"representatives={len(insert_rep_ids)}; clusters={len(insert_summary)}",
    )

    orf_nt_ids = fasta_ids(args.orf_fna)
    orf_aa_ids = fasta_ids(args.orf_faa)
    check(
        "orf_nucleotide_protein_id_consistency",
        orf_nt_ids == orf_aa_ids
        and len(orf_nt_ids) == len(set(orf_nt_ids))
        and (args.expected_orfs == 0 or len(orf_nt_ids) == args.expected_orfs),
        f"nucleotide={len(orf_nt_ids)}; protein={len(orf_aa_ids)}; expected={args.expected_orfs}",
    )

    mapping = {row["orf_id"]: row for row in table_rows(args.orf_mapping)}
    missing_rows = table_rows(args.missing_mappings)
    mapping_complete = set(mapping) == set(orf_nt_ids) and not missing_rows
    mapping_complete &= all(
        row["parent_insert_rep"]
        and row["insert_cluster_id"]
        and row["insert_cluster_size"]
        and row["insert_cluster_id"] in insert_summary
        and row["parent_insert_rep"] in insert_summary
        and int(row["insert_cluster_size"]) == int(insert_summary[row["insert_cluster_id"]]["size"])
        for row in mapping.values()
    )
    check(
        "orf_to_parent_insert_mapping_complete",
        mapping_complete,
        f"orfs={len(orf_nt_ids)}; mapped={len(mapping)}; missing_rows={len(missing_rows)}",
    )

    orf_pairs = cluster_pairs(args.orf_members)
    orf_member_counts = Counter(member for _, member in orf_pairs)
    orf_clusters: dict[str, list[str]] = defaultdict(list)
    for representative, member in orf_pairs:
        orf_clusters[representative].append(member)
    check(
        "orf_cluster_membership_complete",
        set(orf_member_counts) == set(orf_nt_ids)
        and all(count == 1 for count in orf_member_counts.values()),
        f"orfs={len(orf_nt_ids)}; unique_members={len(orf_member_counts)}; rows={len(orf_pairs)}",
    )
    orf_summary_rows = table_rows(args.orf_summary)
    orf_summary = {row["cluster_id"]: row for row in orf_summary_rows}
    orf_summary_consistent = set(orf_summary) == set(orf_clusters)
    if orf_summary_consistent:
        for cluster, row in orf_summary.items():
            members = [member for member in row["members"].split(",") if member]
            orf_summary_consistent &= int(row["size"]) == len(members)
            orf_summary_consistent &= set(members) == set(orf_clusters[cluster])
    check(
        "orf_cluster_summary_consistent",
        orf_summary_consistent,
        f"membership_clusters={len(orf_clusters)}; summary_clusters={len(orf_summary)}",
    )

    support_rows = table_rows(args.support_table)
    support_by_id = {row["reporf_id"]: row for row in support_rows}
    representative_orf_ids = fasta_ids(args.representative_orfs)
    check(
        "reporf_count",
        len(support_rows) == len(orf_clusters) and (args.expected_reporfs == 0 or len(support_rows) == args.expected_reporfs),
        f"support_rows={len(support_rows)}; clusters={len(orf_clusters)}; expected={args.expected_reporfs}",
    )
    check(
        "representative_orf_fasta_complete",
        set(representative_orf_ids) == set(support_by_id)
        and len(representative_orf_ids) == len(support_rows),
        f"representative_fasta={len(representative_orf_ids)}; support_rows={len(support_rows)}",
    )

    support_errors = 0
    for representative, members in orf_clusters.items():
        row = support_by_id.get(representative)
        if row is None:
            support_errors += 1
            continue
        parent_ids = {mapping[member]["parent_insert_rep"] for member in members}
        insert_cluster_ids = {mapping[member]["insert_cluster_id"] for member in members}
        expected_support = sum(int(insert_summary[cluster]["size"]) for cluster in insert_cluster_ids)
        if (
            int(row["orf_cluster_member_count"]) != len(members)
            or int(row["unique_parent_insert_count"]) != len(parent_ids)
            or int(row["unique_insert_cluster_count"]) != len(insert_cluster_ids)
            or int(row["reporf_support"]) != expected_support
            or set(filter(None, row["insert_cluster_ids"].split(","))) != insert_cluster_ids
            or set(filter(None, row["member_orfs"].split(","))) != set(members)
        ):
            support_errors += 1
    check(
        "reporf_support_recalculation",
        support_errors == 0,
        f"support_mismatches={support_errors}; unique_parent_clusters_counted_once=true",
    )

    ranking_keys = [
        (-int(row["reporf_support"]), -int(row["orf_cluster_member_count"]), row["reporf_id"])
        for row in support_rows
    ]
    sequential_ranks = [int(row["reporf_rank"]) for row in support_rows] == list(
        range(1, len(support_rows) + 1)
    )
    check(
        "deterministic_reporf_ranking",
        ranking_keys == sorted(ranking_keys) and sequential_ranks,
        "ranking=support_desc,member_count_desc,reporf_id_asc; ranks_are_sequential="
        + str(sequential_ranks).lower(),
    )

    metrics = selection_metrics(args.selection_summary)
    nominal_n = int(metrics["nominal_n"])
    boundary_support = int(metrics["nominal_boundary_support"])
    primary_rows = table_rows(args.primary_table)
    exploratory_rows = table_rows(args.exploratory_table)
    if metrics.get("selection_mode") == "minimum_support":
        expected_primary_ids = [row["reporf_id"] for row in support_rows if int(row["reporf_support"]) >= int(metrics["minimum_primary_support"])]
    elif boundary_support > 1:
        expected_primary_ids = [row["reporf_id"] for row in support_rows[:nominal_n]]
    else:
        minimum_support = int(metrics["minimum_primary_support"])
        expected_primary_ids = [
            row["reporf_id"]
            for row in support_rows
            if int(row["reporf_support"]) >= minimum_support
        ]
    observed_primary_ids = [row["reporf_id"] for row in primary_rows]
    expected_exploratory_ids = [
        row["reporf_id"] for row in support_rows if int(row["reporf_support"]) == 1
    ]
    observed_exploratory_ids = [row["reporf_id"] for row in exploratory_rows]
    check(
        "tie_aware_selection_rule",
        observed_primary_ids == expected_primary_ids
        and observed_exploratory_ids == expected_exploratory_ids,
        f"nominal_n={nominal_n}; boundary_support={boundary_support}; "
        f"strategy={metrics['selection_strategy']}",
    )
    check(
        "primary_candidate_count",
        (args.expected_primary == 0 or len(primary_rows) == args.expected_primary),
        f"observed={len(primary_rows)}; expected={args.expected_primary}",
    )

    primary_fna_ids = fasta_ids(args.primary_fna)
    primary_faa_ids = fasta_ids(args.primary_faa)
    exploratory_fna_ids = fasta_ids(args.exploratory_fna)
    exploratory_faa_ids = fasta_ids(args.exploratory_faa)
    check(
        "primary_candidate_table_fasta_consistency",
        set(observed_primary_ids) == set(primary_fna_ids) == set(primary_faa_ids)
        and len(observed_primary_ids) == len(set(observed_primary_ids))
        and len(primary_fna_ids) == len(set(primary_fna_ids))
        and len(primary_faa_ids) == len(set(primary_faa_ids)),
        f"table={len(observed_primary_ids)}; fna={len(primary_fna_ids)}; faa={len(primary_faa_ids)}",
    )
    check(
        "exploratory_support1_table_fasta_consistency",
        set(observed_exploratory_ids) == set(exploratory_fna_ids) == set(exploratory_faa_ids)
        and len(observed_exploratory_ids) == len(set(observed_exploratory_ids))
        and len(exploratory_fna_ids) == len(set(exploratory_fna_ids))
        and len(exploratory_faa_ids) == len(set(exploratory_faa_ids)),
        f"table={len(observed_exploratory_ids)}; fna={len(exploratory_fna_ids)}; "
        f"faa={len(exploratory_faa_ids)}",
    )

    annotation_rows = table_rows(args.annotations)
    annotation_ids = [row["reporf_id"] for row in annotation_rows]
    check(
        "all_primary_candidates_have_annotation_rows",
        annotation_ids == observed_primary_ids and len(annotation_ids) == len(set(annotation_ids)),
        f"primary={len(observed_primary_ids)}; annotation_rows={len(annotation_ids)}; "
        f"unique_annotation_ids={len(set(annotation_ids))}",
    )
    retained_no_hits = sum(
        int(row["rgi_hit_count"]) == 0 and int(row["swissprot_hit_count"]) == 0
        for row in annotation_rows
    )
    retained_partial = sum(row["partial_flag"] not in ("", "00") for row in annotation_rows)
    check(
        "annotation_retention_policy",
        len(annotation_rows) == len(primary_rows),
        f"retained_no_hit_candidates={retained_no_hits}; retained_partial_orfs={retained_partial}",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "status", "detail"])
        writer.writerows(checks)
    failures = [row for row in checks if row[1] == "FAIL"]
    print(f"Validation checks={len(checks)}; failures={len(failures)}")
    if failures:
        for name, _, detail in failures:
            print(f"FAILED {name}: {detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
