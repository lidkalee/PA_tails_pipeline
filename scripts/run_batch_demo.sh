#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Prepare the tiny references and verify the legacy one-sample route too.
bash scripts/run_demo.sh
mkdir -p demo/reads/R1 demo/reads/R2
for sample in demoA demoB; do
  cp demo/fixtures/demo_R1.fq "demo/reads/R1/${sample}_R1.fq"
  cp demo/fixtures/demo_R2.fq "demo/reads/R2/${sample}_R2.fq"
done
python scripts/make_samplesheet.py \
  --r1-dir demo/reads/R1 --r2-dir demo/reads/R2 \
  --output demo/batch_samples.tsv --force
snakemake --configfile demo/batch_config.yaml --cores 2 "$@"
Rscript -e 'for (s in c("demoA", "demoB")) {x <- readRDS(paste0("R_output/", s, "_final_big.rds")); stopifnot(nrow(x) == 5, all(x$tail_len == 20))}; cat("Batch demo passed: two independent samples\n")'
