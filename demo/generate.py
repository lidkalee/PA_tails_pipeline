"""Generate a tiny deterministic paired-end dataset for a full workflow smoke test."""

from pathlib import Path
import random
import subprocess

root = Path(__file__).resolve().parent
refs = root / "reference"
fixtures = root / "fixtures"
refs.mkdir(exist_ok=True)
fixtures.mkdir(exist_ok=True)
rng = random.Random(47)
genome = "".join(rng.choices("ACGT", k=5000))
transcript = genome[500:650]
# The tail has 20 T bases; genomic sequence then maps with a 5' soft clip.
r1 = transcript[:100]
for name in ("R1", "R2"):
    with open(fixtures / f"demo_{name}.fq", "w") as out:
        for i in range(5):
            umi = f"{i:015b}".replace("0", "A").replace("1", "C")
            r2 = umi + "GTCAG" + "T" * 20 + genome[1000:1090]
            seq = r1 if name == "R1" else r2
            out.write(f"@demo_read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")
annotation = refs / "annotation.rds"
subprocess.run([
    "Rscript", "-e",
    "saveRDS(data.frame(gene_id='ENSGDEMO', seqnames='chrDemo', "
    "start=950L, end=1150L, strand='+'), commandArgs(TRUE)[1])",
    str(annotation),
], check=True)
print("Generated demo fixtures and annotation")
