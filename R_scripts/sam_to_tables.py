"""Extract columns consumed by the R analysis from SAM files."""

import argparse
import re


def normalized_name(name):
    return re.sub(r"/[12](?=:[ACGTN]{15}$|$)", "", name)


def rows(path):
    with open(path) as handle:
        for line in handle:
            if not line.startswith("@"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(f"Malformed SAM alignment in {path}")
                yield fields


def convert(r1, r2, out_r1, out_r2):
    with open(out_r1, "w") as out:
        for f in rows(r1):
            parts = f[2].split("|")
            if len(parts) < 2:
                raise ValueError("Transcriptome FASTA headers must contain transcript_id|gene_id")
            out.write("\t".join((normalized_name(f[0]), parts[1], parts[0], f[1], f[3], f[4], f[5])) + "\n")
    with open(out_r2, "w") as out:
        for f in rows(r2):
            out.write("\t".join((normalized_name(f[0]), f[1], f[2], f[3], f[5], f[9])) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("r1", "r2", "out-r1", "out-r2"):
        parser.add_argument("--" + key, required=True)
    args = parser.parse_args()
    convert(args.r1, args.r2, args.out_r1, args.out_r2)
