"""Catch STAR builds that silently ignore FASTQ input on Apple Silicon."""

import sys

fastq, sam = sys.argv[1:]
with open(fastq) as handle:
    has_reads = bool(handle.read(1))
if has_reads:
    with open(sam) as handle:
        if not any(not line.startswith("@") for line in handle):
            raise SystemExit(
                "STAR returned no SAM records for nonempty FASTQ. On Apple Silicon, "
                "run bash scripts/install_star_macos.sh; otherwise inspect the STAR log."
            )
