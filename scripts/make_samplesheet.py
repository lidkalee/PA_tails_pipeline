"""Pair FASTQ files from separate R1/R2 folders by a shared sample ID."""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.sample_sheet import SAFE_SAMPLE

DEFAULT_R1 = r"^(?P<sample>.+)[_.-]R1(?:[_-]\d+)?\.(?:fastq|fq)(?:\.gz)?$"
DEFAULT_R2 = r"^(?P<sample>.+)[_.-]R2(?:[_-]\d+)?\.(?:fastq|fq)(?:\.gz)?$"
FASTQ = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)


def _scan(folder, expression, mate):
    folder = Path(folder)
    if not folder.is_dir():
        raise ValueError(f"{mate} directory does not exist: {folder}")
    pattern = re.compile(expression, re.IGNORECASE)
    if "sample" not in pattern.groupindex:
        raise ValueError(f"{mate} regex must have a named (?P<sample>...) group")
    found = {}
    unmatched = []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or not FASTQ.search(path.name):
            continue
        match = pattern.fullmatch(path.name)
        if not match:
            unmatched.append(path.name)
            continue
        sample = match.group("sample")
        if not SAFE_SAMPLE.fullmatch(sample):
            raise ValueError(f"Unsafe sample ID {sample!r} from {path.name}")
        if sample in found:
            raise ValueError(f"Multiple {mate} FASTQs for {sample!r}: {found[sample]} and {path}")
        found[sample] = path.resolve()
    if unmatched:
        raise ValueError(
            f"{mate} FASTQs did not match the regex: {', '.join(unmatched[:5])}. "
            "Use --r1-regex/--r2-regex or write the sample sheet manually."
        )
    if not found:
        raise ValueError(f"No {mate} FASTQs found in {folder}")
    return found


def discover_pairs(r1_dir, r2_dir, r1_regex=DEFAULT_R1, r2_regex=DEFAULT_R2):
    r1 = _scan(r1_dir, r1_regex, "R1")
    r2 = _scan(r2_dir, r2_regex, "R2")
    if r1.keys() != r2.keys():
        only_r1 = ", ".join(sorted(r1.keys() - r2.keys())) or "none"
        only_r2 = ", ".join(sorted(r2.keys() - r1.keys())) or "none"
        raise ValueError(f"Unpaired samples (R1 only: {only_r1}; R2 only: {only_r2})")
    return {sample: (r1[sample], r2[sample]) for sample in sorted(r1)}


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
    parser.add_argument("--r1-dir", required=True, help="Directory containing R1 FASTQs")
    parser.add_argument("--r2-dir", required=True, help="Directory containing R2 FASTQs")
    parser.add_argument("--output", default="samples.tsv", help="Output TSV (default: samples.tsv)")
    parser.add_argument("--r1-regex", default=DEFAULT_R1, help="Regex with named sample group")
    parser.add_argument("--r2-regex", default=DEFAULT_R2, help="Regex with named sample group")
    parser.add_argument("--force", action="store_true", help="Replace an existing output sheet")
    args = parser.parse_args()
    try:
        pairs = discover_pairs(args.r1_dir, args.r2_dir, args.r1_regex, args.r2_regex)
        write_sheet(pairs, args.output, args.force)
    except (ValueError, FileExistsError, re.error) as error:
        parser.error(str(error))
    print(f"Wrote {len(pairs)} paired samples to {args.output}")


if __name__ == "__main__":
    main()
