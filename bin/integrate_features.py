#!/usr/bin/env python3
"""Join Prodigal, ABRicate, and ISEScan records by contig and coordinates."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


FEATURE_FIELDS = [
    "sample",
    "contig",
    "feature_type",
    "start",
    "end",
    "length_bp",
    "strand",
    "label",
    "identity_pct",
    "reference_coverage_pct",
    "database",
    "accession",
    "product",
    "completeness",
    "attributes",
]
PAIR_FIELDS = [
    "sample",
    "contig",
    "feature_1_type",
    "feature_1_label",
    "feature_1_start",
    "feature_1_end",
    "feature_1_strand",
    "feature_2_type",
    "feature_2_label",
    "feature_2_start",
    "feature_2_end",
    "feature_2_strand",
    "gap_bp",
    "overlap_bp",
    "orientation",
    "pair_class",
]
CONTEXT_FIELDS = [
    "sample",
    "contig",
    "hit_type",
    "hit_label",
    "hit_start",
    "hit_end",
    "hit_identity_pct",
    "hit_reference_coverage_pct",
    "cds_id",
    "cds_start",
    "cds_end",
    "cds_strand",
    "cds_product",
    "distance_bp",
    "overlaps_hit",
]


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    length = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"): 
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def parse_prodigal(path: Path, sample: str) -> list[dict]:
    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = parse_attributes(cols[8])
            start, end = int(cols[3]), int(cols[4])
            records.append(
                {
                    "sample": sample,
                    "contig": cols[0],
                    "feature_type": "CDS",
                    "start": start,
                    "end": end,
                    "length_bp": abs(end - start) + 1,
                    "strand": cols[6],
                    "label": attrs.get("ID", attrs.get("gene", "Prodigal_CDS")),
                    "identity_pct": "",
                    "reference_coverage_pct": "",
                    "database": "Prodigal",
                    "accession": "",
                    "product": attrs.get("product", ""),
                    "completeness": "",
                    "attributes": cols[8],
                }
            )
    return records


def get(row: dict[str, str], *keys: str, default: str = "") -> str:
    lowered = {str(k).lstrip("#").lower(): v for k, v in row.items() if k is not None}
    for key in keys:
        value = lowered.get(key.lower())
        if value is not None:
            return value.strip()
    return default


def numeric(value: str, default: float = 0.0) -> float:
    try:
        return float(value.rstrip("%"))
    except (TypeError, ValueError):
        return default


def parse_abricate(path: Path, sample: str, database: str) -> list[dict]:
    records = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig = get(row, "SEQUENCE")
            start = int(get(row, "START", default="0") or 0)
            end = int(get(row, "END", default="0") or 0)
            if not contig or start < 1 or end < 1:
                continue
            db = get(row, "DATABASE", default=database).lower() or database
            feature_type = "AMR" if db == "resfinder" else "VFDB"
            records.append(
                {
                    "sample": sample,
                    "contig": contig,
                    "feature_type": feature_type,
                    "start": min(start, end),
                    "end": max(start, end),
                    "length_bp": abs(end - start) + 1,
                    "strand": get(row, "STRAND"),
                    "label": get(row, "GENE", default="unlabelled"),
                    "identity_pct": get(row, "%IDENTITY"),
                    "reference_coverage_pct": get(row, "%COVERAGE"),
                    "database": db,
                    "accession": get(row, "ACCESSION"),
                    "product": get(row, "PRODUCT"),
                    "completeness": "",
                    "attributes": "resistance=" + get(row, "RESISTANCE"),
                }
            )
    return records


def parse_isescan(path: Path, sample: str) -> list[dict]:
    records = []
    if not path.exists() or path.stat().st_size == 0:
        return records
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig = get(row, "seqID", "seqid", "sequence")
            start_text = get(row, "isBegin", "start", "begin")
            end_text = get(row, "isEnd", "end", "stop")
            if not contig or not start_text or not end_text:
                continue
            start, end = int(float(start_text)), int(float(end_text))
            completeness = get(row, "type", "completeness")
            if completeness == "c":
                completeness = "complete"
            elif completeness == "p":
                completeness = "partial"
            records.append(
                {
                    "sample": sample,
                    "contig": contig,
                    "feature_type": "IS",
                    "start": min(start, end),
                    "end": max(start, end),
                    "length_bp": abs(end - start) + 1,
                    "strand": get(row, "strand"),
                    "label": get(row, "family", "cluster", default="ISEScan_call"),
                    "identity_pct": "",
                    "reference_coverage_pct": "",
                    "database": "ISEScan",
                    "accession": get(row, "cluster"),
                    "product": "Insertion sequence candidate",
                    "completeness": completeness,
                    "attributes": "\t".join(f"{k}={v}" for k, v in row.items() if v not in (None, "")),
                }
            )
    return records


def interval_gap(a: dict, b: dict) -> tuple[int, int]:
    overlap = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    gap = max(0, max(a["start"], b["start"]) - min(a["end"], b["end"]) - 1)
    return gap, overlap


def orientation(a: dict, b: dict) -> str:
    if a["strand"] not in ("+", "-") or b["strand"] not in ("+", "-"):
        return "unknown"
    return "same" if a["strand"] == b["strand"] else "opposite"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--prodigal-gff", required=True, type=Path)
    parser.add_argument("--resfinder", required=True, type=Path)
    parser.add_argument("--vfdb", required=True, type=Path)
    parser.add_argument("--isescan-tsv", required=True, type=Path)
    parser.add_argument("--isescan-gff", required=True, type=Path)
    parser.add_argument("--max-distance", required=True, type=int)
    parser.add_argument("--min-identity", required=True, type=float)
    parser.add_argument("--min-coverage", required=True, type=float)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    contig_lengths = read_fasta_lengths(args.assembly)
    features = []
    features.extend(parse_prodigal(args.prodigal_gff, args.sample))
    features.extend(parse_abricate(args.resfinder, args.sample, "resfinder"))
    features.extend(parse_abricate(args.vfdb, args.sample, "vfdb"))
    features.extend(parse_isescan(args.isescan_tsv, args.sample))
    features.sort(key=lambda x: (x["contig"], x["start"], x["end"], x["feature_type"], x["label"]))
    write_tsv(Path(f"{args.prefix}.contig_features.tsv"), FEATURE_FIELDS, features)

    candidates = [f for f in features if f["feature_type"] in {"AMR", "VFDB", "IS"}]
    by_contig: dict[str, list[dict]] = defaultdict(list)
    genes_by_contig: dict[str, list[dict]] = defaultdict(list)
    hits = []
    for feature in features:
        if feature["feature_type"] == "CDS":
            genes_by_contig[feature["contig"]].append(feature)
        elif feature["feature_type"] in {"AMR", "VFDB"}:
            hits.append(feature)
        if feature["feature_type"] in {"AMR", "VFDB", "IS"}:
            by_contig[feature["contig"]].append(feature)

    pairs = []
    for contig, items in by_contig.items():
        items.sort(key=lambda x: (x["start"], x["end"], x["feature_type"], x["label"]))
        for index, first in enumerate(items):
            for second in items[index + 1 :]:
                if second["start"] - first["end"] - 1 > args.max_distance:
                    break
                gap, overlap = interval_gap(first, second)
                if gap > args.max_distance:
                    continue
                left, right = (first, second) if (first["start"], first["end"]) <= (second["start"], second["end"]) else (second, first)
                pair_types = sorted([first["feature_type"], second["feature_type"]])
                pairs.append(
                    {
                        "sample": args.sample,
                        "contig": contig,
                        "feature_1_type": left["feature_type"],
                        "feature_1_label": left["label"],
                        "feature_1_start": left["start"],
                        "feature_1_end": left["end"],
                        "feature_1_strand": left["strand"],
                        "feature_2_type": right["feature_type"],
                        "feature_2_label": right["label"],
                        "feature_2_start": right["start"],
                        "feature_2_end": right["end"],
                        "feature_2_strand": right["strand"],
                        "gap_bp": gap,
                        "overlap_bp": overlap,
                        "orientation": orientation(left, right),
                        "pair_class": "-".join(pair_types),
                    }
                )
    pairs.sort(key=lambda x: (x["contig"], x["feature_1_start"], x["feature_2_start"]))
    write_tsv(Path(f"{args.prefix}.nearby_locus_pairs.tsv"), PAIR_FIELDS, pairs)

    contexts = []
    for hit in hits:
        nearby = []
        for gene in genes_by_contig.get(hit["contig"], []):
            gap, overlap = interval_gap(hit, gene)
            if gap <= args.max_distance:
                nearby.append((gap, -overlap, gene))
        for gap, neg_overlap, gene in sorted(nearby, key=lambda x: (x[0], x[1], x[2]["start"])):
            contexts.append(
                {
                    "sample": args.sample,
                    "contig": hit["contig"],
                    "hit_type": hit["feature_type"],
                    "hit_label": hit["label"],
                    "hit_start": hit["start"],
                    "hit_end": hit["end"],
                    "hit_identity_pct": hit["identity_pct"],
                    "hit_reference_coverage_pct": hit["reference_coverage_pct"],
                    "cds_id": gene["label"],
                    "cds_start": gene["start"],
                    "cds_end": gene["end"],
                    "cds_strand": gene["strand"],
                    "cds_product": gene["product"],
                    "distance_bp": gap,
                    "overlaps_hit": str(neg_overlap < 0).lower(),
                }
            )
    write_tsv(Path(f"{args.prefix}.hit_gene_context.tsv"), CONTEXT_FIELDS, contexts)

    counts = Counter(feature["feature_type"] for feature in features)
    summary = {
        "sample": args.sample,
        "assembly": args.assembly.name,
        "contigs": len(contig_lengths),
        "assembly_bp": sum(contig_lengths.values()),
        "prodigal_cds": counts.get("CDS", 0),
        "abricate_resfinder_hits": counts.get("AMR", 0),
        "abricate_vfdb_hits": counts.get("VFDB", 0),
        "isescan_elements": counts.get("IS", 0),
        "candidate_feature_pairs_within_max_distance": len(pairs),
        "hit_gene_context_rows": len(contexts),
        "max_pair_distance_bp": args.max_distance,
        "abricate_min_identity_pct": args.min_identity,
        "abricate_min_reference_coverage_pct": args.min_coverage,
        "isescan_gff_available": args.isescan_gff.exists() and args.isescan_gff.stat().st_size > 0,
    }
    Path(f"{args.prefix}.summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
