"""Label Bowtie2 genome hits unique or multi using up to two reported alignments."""

import argparse
from collections import Counter


def parse_sam(sam_path, output_path):
    counts = Counter()
    with open(sam_path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed SAM alignment in {sam_path}")
            if not (int(fields[1]) & 4):
                counts[fields[0]] += 1
    with open(output_path, "w") as handle:
        for name, count in counts.items():
            handle.write(f"{name}\t{'unique' if count == 1 else 'multi'}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    parse_sam(args.input, args.output)
