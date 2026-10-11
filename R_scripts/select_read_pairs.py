"""Keep synchronized pairs whose R1 aligned and whose R2 has the adapter motif."""

import argparse
import gzip
import re
from itertools import zip_longest

ADAPTER = re.compile(r"^[ACGTN]{15}GTCAG", re.IGNORECASE)


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def records(handle, path):
    while True:
        lines = [handle.readline() for _ in range(4)]
        if not lines[0]:
            return
        if any(not line for line in lines) or not lines[0].startswith("@") or not lines[2].startswith("+"):
            raise ValueError(f"Malformed FASTQ record in {path}")
        if len(lines[1].rstrip("\r\n")) != len(lines[3].rstrip("\r\n")):
            raise ValueError(f"Sequence and quality lengths differ in {path}")
        yield lines


def read_name(header):
    name = header[1:].split()[0]
    return re.sub(r"/[12]$", "", name)


def select(sam, r1, r2, out_r1, out_r2):
    aligned = set()
    with open(sam) as handle:
        for line in handle:
            if not line.startswith("@"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and not (int(fields[1]) & 4):
                    aligned.add(read_name("@" + fields[0]))
    count = 0
    with open_fastq(r1) as h1, open_fastq(r2) as h2, open(out_r1, "w") as o1, open(out_r2, "w") as o2:
        for a, b in zip_longest(records(h1, r1), records(h2, r2)):
            if a is None or b is None:
                raise ValueError("R1 and R2 contain different numbers of records")
            name = read_name(a[0])
            if name != read_name(b[0]):
                raise ValueError(f"R1 and R2 are out of sync at {name}")
            if name in aligned and ADAPTER.match(b[1].strip()):
                o1.writelines(a)
                o2.writelines(b)
                count += 1
    print(f"Selected {count} read pairs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("sam", "r1", "r2", "out-r1", "out-r2"):
        parser.add_argument("--" + key, required=True)
    args = parser.parse_args()
    select(args.sam, args.r1, args.r2, args.out_r1, args.out_r2)
