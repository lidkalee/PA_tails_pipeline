#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for command in python Rscript bowtie2-build fastp snakemake; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "Missing $command. Create and activate environment.yml first." >&2
    exit 1
  fi
done
if [[ "$(uname -s)" == Darwin && "$(uname -m)" == arm64 && ! -x tools/STAR ]]; then
  bash scripts/install_star_macos.sh
fi
python demo/generate.py
if [[ ! -f demo/reference/transcript_index.1.bt2 && ! -f demo/reference/transcript_index.1.bt2l ]]; then
  bowtie2-build --threads 2 demo/fixtures/transcripts.fa demo/reference/transcript_index >/dev/null
fi
if [[ ! -f demo/reference/genome_index.1.bt2 && ! -f demo/reference/genome_index.1.bt2l ]]; then
  bowtie2-build --threads 2 demo/fixtures/genome.fa demo/reference/genome_index >/dev/null
fi
if [[ ! -f demo/reference/star_index/Genome ]]; then
  mkdir -p demo/reference/star_index
  bash scripts/star.sh --runMode genomeGenerate --runThreadN 2 \
    --genomeDir demo/reference/star_index \
    --genomeFastaFiles demo/fixtures/genome.fa \
    --genomeSAindexNbases 3 --genomeChrBinNbits 8
fi
snakemake --configfile demo/config.yaml --cores 2 "$@"
Rscript -e 'x <- readRDS("R_output/demo_final_big.rds"); d <- read.csv("R_output/demo_deduplication.csv"); stopifnot(nrow(x) == 5, all(x$geneR1 == "ENSGDEMO"), all(x$tail_len == 20), all(x$type == "cigar"), length(unique(x$UMI)) == 5, identical(d$reads, c(5L, 5L))); cat("Demo passed: five annotated reads, each with a 20-base tail\n")'
