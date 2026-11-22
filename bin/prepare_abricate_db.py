#!/usr/bin/env python3
"""Convert a pinned nucleotide reference FASTA to ABRicate custom DB format."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def fasta_records(path: Path):
    header = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif header is None:
                raise ValueError(f"Sequence appears before FASTA header in {path}")
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def header_fields(database: str, original: str) -> tuple[str, str, str]:
    first, _, description = original.partition(" ")
    description = description.strip() or original
    if database == "vfdb":
        # VFDB Set A identifiers often look like VFG... (gb|WP_...) (gene) ...
        gene = first.split("(", 1)[0]
        accession_match = re.search(r"\((?:gb\|)?([^()|]+)(?:\|[^()]*)?\)", first)
        accession = accession_match.group(1) if accession_match else first
    else:
        # Preserve the full reference token: ResFinder headers commonly combine
        # an allele name and accession in a single token.
        gene = first
        accession = first
    safe = lambda value: re.sub(r"~~~", "_", value).replace("\t", " ").strip()
    return safe(gene), safe(accession), safe(description)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True, choices=("resfinder", "vfdb"))
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    args = parser.parse_args()

    digest = hashlib.sha256()
    count = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out, args.metadata.open(
        "w", encoding="utf-8", newline=""
    ) as meta:
        writer = csv.writer(meta, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["database", "record_number", "original_header", "abricate_gene", "accession", "description"]
        )
        for number, (header, sequence) in enumerate(fasta_records(args.input), start=1):
            if not sequence:
                continue
            gene, accession, description = header_fields(args.database, header)
            # ABRicate custom databases encode the DB, gene, accession, and
            # resistance/product fields in the FASTA description.
            product = re.sub(r"~~~", "_", description).replace("\t", " ").strip()
            out.write(f">{args.database}~~~{gene}~~~{accession}~~~NA {product}\n")
            for offset in range(0, len(sequence), 80):
                line = sequence[offset : offset + 80]
                out.write(line + "\n")
                digest.update(line.upper().encode("ascii"))
            count += 1
            writer.writerow([args.database, number, header, gene, accession, description])

    if count == 0:
        raise ValueError(f"No FASTA records found in {args.input}")
    with args.metadata.open("a", encoding="utf-8") as meta:
        meta.write(f"# records={count}\n")
        meta.write(f"# normalized_sequence_sha256={digest.hexdigest()}\n")
        source_digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
        meta.write(f"# source_fasta_sha256={source_digest}\n")


if __name__ == "__main__":
    main()
