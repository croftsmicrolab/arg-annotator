#!/usr/bin/env python3
"""Merge support, ORF completeness, CARD/RGI, and Swiss-Prot evidence without filtering."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path


DIAMOND_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qcovhsp",
    "scovhsp",
    "stitle",
]


def number(value, default=math.nan):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def first(row, *names):
    for name in names:
        value = row.get(name, "")
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def clean_orf(value: str) -> str:
    cleaned = (value or "").split("#", 1)[0].strip().split()
    return cleaned[0] if cleaned else ""


def fasta_lengths(path: Path):
    lengths = {}
    current = None
    length = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = clean_orf(line[1:])
                length = 0
            else:
                length += len(line.strip().rstrip("*"))
        if current is not None:
            lengths[current] = length
    return lengths


def parse_swissprot_title(title: str, subject: str):
    accession = subject.split("|")[1] if subject.count("|") >= 2 else subject.split()[0]
    cleaned = re.sub(r"^(?:sp|tr)\|[^|]+\|\S+\s+", "", title.strip())
    organism_match = re.search(r"\bOS=(.+?)(?:\sOX=|\sGN=|\sPE=|\sSV=|$)", cleaned)
    organism = organism_match.group(1).strip() if organism_match else ""
    protein = re.split(r"\sOS=", cleaned, maxsplit=1)[0].strip()
    return accession, protein, organism


def read_diamond(path: Path):
    best = {}
    hit_counts = Counter()
    if not path.exists() or path.stat().st_size == 0:
        return best, hit_counts
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=DIAMOND_FIELDS)
        for row in reader:
            query = clean_orf(row["qseqid"])
            hit_counts[query] += 1
            key = (
                number(row["evalue"], math.inf),
                -number(row["bitscore"], -math.inf),
                -number(row["pident"], -math.inf),
                row["sseqid"],
            )
            if query in best and key >= best[query][0]:
                continue
            accession, protein, organism = parse_swissprot_title(row["stitle"], row["sseqid"])
            best[query] = (
                key,
                {
                    "swissprot_accession": accession,
                    "swissprot_protein_name": protein,
                    "swissprot_organism": organism,
                    "swissprot_percent_identity": row["pident"],
                    "swissprot_alignment_length": row["length"],
                    "swissprot_query_coverage": row["qcovhsp"],
                    "swissprot_subject_coverage": row["scovhsp"],
                    "swissprot_evalue": row["evalue"],
                    "swissprot_bitscore": row["bitscore"],
                    "swissprot_description": row["stitle"],
                },
            )
    return {query: value[1] for query, value in best.items()}, hit_counts


def cutoff_priority(value: str) -> int:
    text = value.lower()
    if "perfect" in text:
        return 0
    if "strict" in text:
        return 1
    if "loose" in text:
        return 2
    return 9


def read_rgi(path: Path, candidate_ids: set[str]):
    best = {}
    hit_counts = Counter()
    if not path.exists() or path.stat().st_size == 0:
        return best, hit_counts
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query = clean_orf(first(row, "Contig", "contig", "sequence_id"))
            if query not in candidate_ids:
                query = clean_orf(first(row, "ORF_ID", "orf_id"))
                while query not in candidate_ids and re.search(r"_[0-9]+$", query):
                    query = re.sub(r"_[0-9]+$", "", query)
            if query not in candidate_ids:
                continue
            hit_counts[query] += 1
            cutoff = first(row, "Cut_Off", "Cut Off", "hit_type")
            key = (
                cutoff_priority(cutoff),
                -number(first(row, "Best_Hit_Bitscore", "Best Hit Bitscore", "bitscore"), -math.inf),
                -number(first(row, "Best_Identities", "Best Identities", "percent_identity"), -math.inf),
                first(row, "Best_Hit_ARO", "Best Hit ARO", "card_hit_name"),
            )
            if query in best and key >= best[query][0]:
                continue
            aro = first(row, "ARO", "ARO Accession", "aro_accession")
            if aro.isdigit():
                aro = f"ARO:{aro}"
            best[query] = (
                key,
                {
                    "rgi_hit_name": first(row, "Best_Hit_ARO", "Best Hit ARO", "card_hit_name"),
                    "rgi_aro_accession": aro,
                    "rgi_hit_type": cutoff,
                    "rgi_model_type": first(row, "Model_type", "Model Type", "model_type"),
                    "rgi_drug_class": first(row, "Drug Class", "Drug_Class", "drug_class"),
                    "rgi_resistance_mechanism": first(
                        row, "Resistance Mechanism", "Resistance_Mechanism", "resistance_mechanism"
                    ),
                    "rgi_amr_gene_family": first(
                        row, "AMR Gene Family", "AMR_Gene_Family", "amr_gene_family"
                    ),
                    "rgi_percent_identity": first(
                        row, "Best_Identities", "Best Identities", "percent_identity"
                    ),
                    "rgi_model_coverage": first(
                        row,
                        "Percentage Length of Reference Sequence",
                        "Model Coverage",
                        "model_coverage",
                    ),
                    "rgi_bitscore": first(
                        row, "Best_Hit_Bitscore", "Best Hit Bitscore", "bitscore"
                    ),
                    "rgi_antibiotic": first(row, "Antibiotic", "antibiotic"),
                },
            )
    return {query: value[1] for query, value in best.items()}, hit_counts


def classify(row, thresholds):
    hit_type = row.get("rgi_hit_type", "").lower()
    protein = row.get("swissprot_protein_name", "").lower()
    identity = number(row.get("swissprot_percent_identity"))
    coverage = number(row.get("swissprot_query_coverage"))
    hypothetical = any(
        token in protein for token in ("hypothetical", "uncharacterized", "unknown protein")
    )

    if "perfect" in hit_type or "strict" in hit_type:
        evidence, confidence = "known_or_close_CARD_homolog", "high"
    elif "loose" in hit_type:
        evidence, confidence = "distant_CARD_homolog", "medium"
    elif row.get("swissprot_accession") and not hypothetical:
        if identity >= thresholds["medium_identity"] and coverage >= thresholds["medium_coverage"]:
            evidence, confidence = "functional_homolog", "medium"
        elif identity >= thresholds["low_identity"] and coverage >= thresholds["low_coverage"]:
            evidence, confidence = "low_similarity_functional_homolog", "low"
        else:
            evidence, confidence = "very_remote_functional_homolog", "exploratory"
    elif row.get("swissprot_accession"):
        evidence, confidence = "hypothetical_or_uncharacterized_homolog", "exploratory"
    else:
        evidence, confidence = "no_CARD_or_SwissProt_hit", "exploratory"
    retention = "established" if confidence == "high" else "novelty_preserved_candidate"
    return evidence, confidence, retention


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--candidate-faa", required=True, type=Path)
    parser.add_argument("--orf-mapping", required=True, type=Path)
    parser.add_argument("--diamond", required=True, type=Path)
    parser.add_argument("--rgi-status", choices=["searched", "disabled"], default="searched")
    parser.add_argument("--swissprot-status", choices=["searched", "disabled"], default="searched")
    parser.add_argument("--rgi", required=True, type=Path)
    parser.add_argument("--medium-identity", required=True, type=float)
    parser.add_argument("--medium-query-coverage", required=True, type=float)
    parser.add_argument("--low-identity", required=True, type=float)
    parser.add_argument("--low-query-coverage", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    with args.candidate_table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        candidates = list(reader)
        base_fields = list(reader.fieldnames or [])
    candidate_ids = {row["reporf_id"] for row in candidates}

    mapping = {}
    with args.orf_mapping.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["orf_id"] in candidate_ids:
                mapping[row["orf_id"]] = row
    lengths = fasta_lengths(args.candidate_faa)
    diamond, diamond_counts = read_diamond(args.diamond)
    rgi, rgi_counts = read_rgi(args.rgi, candidate_ids)
    thresholds = {
        "medium_identity": args.medium_identity,
        "medium_coverage": args.medium_query_coverage,
        "low_identity": args.low_identity,
        "low_coverage": args.low_query_coverage,
    }

    annotation_fields = [
        "rgi_search_status",
        "swissprot_search_status",
        "protein_length",
        "orf_start",
        "orf_end",
        "orf_strand",
        "partial_flag",
        "rgi_hit_count",
        "rgi_hit_name",
        "rgi_aro_accession",
        "rgi_hit_type",
        "rgi_model_type",
        "rgi_drug_class",
        "rgi_resistance_mechanism",
        "rgi_amr_gene_family",
        "rgi_percent_identity",
        "rgi_model_coverage",
        "rgi_bitscore",
        "rgi_antibiotic",
        "swissprot_hit_count",
        "swissprot_accession",
        "swissprot_protein_name",
        "swissprot_organism",
        "swissprot_percent_identity",
        "swissprot_alignment_length",
        "swissprot_query_coverage",
        "swissprot_subject_coverage",
        "swissprot_evalue",
        "swissprot_bitscore",
        "swissprot_description",
        "evidence_class",
        "annotation_confidence",
        "candidate_retention_class",
    ]
    output_rows = []
    for base in candidates:
        query = base["reporf_id"]
        metadata = mapping.get(query, {})
        row = {
            **base,
            "rgi_search_status": args.rgi_status,
            "swissprot_search_status": args.swissprot_status,
            "protein_length": lengths.get(query, ""),
            "orf_start": metadata.get("start", ""),
            "orf_end": metadata.get("end", ""),
            "orf_strand": metadata.get("strand", ""),
            "partial_flag": metadata.get("partial_flag", ""),
            "rgi_hit_count": rgi_counts.get(query, 0),
            **rgi.get(query, {}),
            "swissprot_hit_count": diamond_counts.get(query, 0),
            **diamond.get(query, {}),
        }
        evidence, confidence, retention = classify(row, thresholds)
        if args.rgi_status == args.swissprot_status == "disabled":
            evidence, confidence, retention = "annotation_not_run", "not_assessed", "retained_candidate"
        elif evidence == "no_CARD_or_SwissProt_hit" and "disabled" in (args.rgi_status, args.swissprot_status):
            evidence = "no_hit_in_searched_database"
        row["evidence_class"] = evidence
        row["annotation_confidence"] = confidence
        row["candidate_retention_class"] = retention
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=base_fields + annotation_fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    metrics = Counter(row["evidence_class"] for row in output_rows)
    metrics.update(
        {
            f"confidence:{key}": value
            for key, value in Counter(row["annotation_confidence"] for row in output_rows).items()
        }
    )
    metrics["total_candidates"] = len(output_rows)
    metrics["partial_orfs"] = sum(
        row.get("partial_flag", "") not in ("", "00") for row in output_rows
    )
    metrics["candidates_with_rgi_hits"] = sum(int(row["rgi_hit_count"]) > 0 for row in output_rows)
    metrics["candidates_with_swissprot_hits"] = sum(
        int(row["swissprot_hit_count"]) > 0 for row in output_rows
    )
    with args.summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["rgi_search_status", args.rgi_status if candidates else "no_candidates"])
        writer.writerow(["swissprot_search_status", args.swissprot_status if candidates else "no_candidates"])
        for key in sorted(metrics):
            writer.writerow([key, metrics[key]])
    print(f"Merged annotations for {len(output_rows)} candidates; no candidates were discarded")


if __name__ == "__main__":
    main()
