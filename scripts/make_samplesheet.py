"""Pair FASTQ files from one or two folders by a shared sample ID."""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.fastq_pairs import DEFAULT_R1, DEFAULT_R2, discover_pairs


def _display_path(path):
    try:
        return str(path.relative_to(Path.cwd()))
    except ValueError:
        return str(path)


def write_sheet(pairs, output, force=False):
    output = Path(output)
    if output.exists() and not force:
        raise FileExistsError(f"{output} already exists; use --force to replace it")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "R1", "R2"))
        for sample, (r1, r2) in pairs.items():
            writer.writerow((sample, _display_path(r1), _display_path(r2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq-dir", help="One directory containing both R1 and R2 FASTQs")
    parser.add_argument("--r1-dir", help="Directory containing R1 FASTQs")
    parser.add_argument("--r2-dir", help="Directory containing R2 FASTQs")
    parser.add_argument("--output", default="samples.tsv", help="Output TSV (default: samples.tsv)")
    parser.add_argument("--r1-regex", default=DEFAULT_R1, help="Regex with named sample group")
    parser.add_argument("--r2-regex", default=DEFAULT_R2, help="Regex with named sample group")
    parser.add_argument("--force", action="store_true", help="Replace an existing output sheet")
    args = parser.parse_args()
    if args.fastq_dir and (args.r1_dir or args.r2_dir):
        parser.error("Use --fastq-dir or both --r1-dir and --r2-dir")
    if args.fastq_dir:
        r1_dir = r2_dir = args.fastq_dir
    elif args.r1_dir and args.r2_dir:
        r1_dir, r2_dir = args.r1_dir, args.r2_dir
    else:
        parser.error("Provide --fastq-dir or both --r1-dir and --r2-dir")
    try:
        pairs = discover_pairs(r1_dir, r2_dir, args.r1_regex, args.r2_regex)
        write_sheet(pairs, args.output, args.force)
    except (ValueError, FileExistsError, re.error) as error:
        parser.error(str(error))
    print(f"Wrote {len(pairs)} paired samples to {args.output}")


if __name__ == "__main__":
    main()
