#!/usr/bin/env python3
"""Copy a possibly gzip-compressed assembly into a plain FASTA file."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    source = gzip.open(args.input, "rb") if compressed else args.input.open("rb")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with source, args.output.open("wb") as target:
        shutil.copyfileobj(source, target, length=4 * 1024 * 1024)


if __name__ == "__main__":
    main()
