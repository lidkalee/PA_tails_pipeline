#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Build the tiny references and verify the legacy one-sample route.
bash scripts/run_demo.sh
mkdir -p demo/reads/R1 demo/reads/R2 demo/reads/mixed
for sample in demoA demoB; do
  cp demo/fixtures/demo_R1.fq "demo/reads/R1/${sample}_R1.fq"
  cp demo/fixtures/demo_R2.fq "demo/reads/R2/${sample}_R2.fq"
done
for sample in demoC demoD; do
  cp demo/fixtures/demo_R1.fq "demo/reads/mixed/${sample}_R1.fq"
  cp demo/fixtures/demo_R2.fq "demo/reads/mixed/${sample}_R2.fq"
done
snakemake --configfile demo/config.yaml \
  --config r1_dir=demo/reads/R1 r2_dir=demo/reads/R2 --cores 2 "$@"
python scripts/make_samplesheet.py \
  --r1-dir demo/reads/R1 --r2-dir demo/reads/R2 \
  --output demo/batch_samples.tsv --force
snakemake --configfile demo/batch_config.yaml --cores 2 "$@"
snakemake --configfile demo/config.yaml \
  --config fastq_dir=demo/reads/mixed --cores 2 "$@"
Rscript -e 'for (s in c("demoA", "demoB", "demoC", "demoD")) {x <- readRDS(paste0("R_output/", s, "_final_big.rds")); stopifnot(nrow(x) == 5, all(x$tail_len == 20))}; cat("Batch demo passed: four independent samples from separate and shared folders\n")'
