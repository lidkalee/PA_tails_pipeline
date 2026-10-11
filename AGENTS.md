# Agent guide

This repository is a Snakemake pipeline for paired-end PA-tail analysis. Keep code, documentation, comments, and user-facing messages in English.

## Working in this repository

- Work from the repository root. Keep input and reference paths configurable in `config.yaml` or through Snakemake `--config` arguments; never commit paths specific to one computer.
- Preserve the legacy single-sample configuration and the batch modes. Batch inputs may be supplied with `fastq_dir`, with `r1_dir` and `r2_dir`, or with a `samples.tsv` file.
- Pair FASTQs by the shared sample name around the `R1`/`R2` token. Reject missing mates, duplicate sample IDs, and ambiguous names rather than guessing.
- Keep generated indexes, intermediate files, and results out of Git. The small files in `demo/fixtures` are test inputs.
- Update `README.md` when installation, configuration, commands, or output paths change.

## Verification

Run these checks after relevant changes:

```bash
python3 -m unittest discover -s tests -v
git diff --check
bash scripts/run_batch_demo.sh
```

The batch demo runs the single-sample example, separate R1/R2 directories, a shared FASTQ directory, and the optional sample-sheet route. It requires the tools from `environment.yml`; see `README.md` for installation and the Apple Silicon STAR setup. GitHub Actions runs the same tests and demo on Linux.
