# PA tails pipeline

This Snakemake workflow analyses A/T-rich tails in paired-end RNA sequencing reads. It aligns R1 to a protein-coding transcriptome, retains matching R2 reads with the expected 15-base UMI followed by `GTCAG`, trims and extracts the UMI with fastp, aligns R2 to an AT-masked genome with STAR, and writes per-read R tables. An additional Bowtie2 genome alignment is used to flag short tails that also align end to end.

The workflow runs on Linux and macOS with a Conda-compatible package manager. Native Windows users can run it under WSL2. On Apple Silicon, the current Bioconda STAR 2.7.11b build can silently read zero records; the repository includes an installer for a pinned upstream fix. The bundled annotation is a human gene annotation RDS; full analysis requires external human reference files and indexes, which are too large for this repository.

## Install

Install [Miniforge](https://github.com/conda-forge/miniforge) or [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), clone this repository, and run from its root:

```bash
conda env create -f environment.yml
conda activate pa-tails
```

With micromamba, use `micromamba create -f environment.yml` and `micromamba activate pa-tails`. The environment contains Snakemake, Bowtie2, STAR, fastp, bedtools, R, and the required R packages. On Apple Silicon, run `bash scripts/install_star_macos.sh` after activation; the demo script does this automatically. This compiles a [pinned upstream STAR fix](https://github.com/alexdobin/STAR/pull/2691) into ignored `tools/STAR` and needs the Conda Clang++ compiler plus Apple Command Line Tools. The workflow uses this local binary on Apple Silicon. The short, version-ranged file can be solved for the current machine architecture; the previous Linux-specific export is intentionally replaced. Bioconda currently provides [Bowtie2](https://bioconda.github.io/recipes/bowtie2/README.html) and [STAR](https://bioconda.github.io/recipes/star/README.html) for Linux and macOS, including ARM variants.

## Run the bundled demo

```bash
bash scripts/run_demo.sh
```

The script uses five committed synthetic read pairs and tiny reference FASTAs in `demo/fixtures`, builds Bowtie2 and STAR indexes, runs the complete workflow, and verifies a detected tail. Generated indexes, annotation, and outputs are ignored by Git. Output appears in `R_output/demo_final_big.rds`, `R_output/demo_final_small.rds`, and `R_output/demo_deduplication.csv`. This is a smoke test of the workflow, not a biological validation dataset.

For fast helper checks without the bioinformatics environment (also run by Linux CI):

```bash
python3 -m unittest discover -s tests -v
```

## Prepare human references

Use matching human reference releases. The original analysis used GENCODE v47 `gencode.v47.pc_transcripts.fa`, `GRCh38.primary_assembly.genome.fa`, and a matching GTF from [GENCODE](https://www.gencodegenes.org/human/). The bundled `annotation/gene_annotations.rds` has `gene_id`, `seqnames`, `start`, `end`, and `strand` columns; gene IDs and chromosome naming in the input references must agree with it. Transcript FASTA headers must begin with `transcript_id|gene_id` (as in the GENCODE protein-coding transcript FASTA). Paths and index prefixes are configurable, so these names are examples rather than required filenames.

Build the transcriptome and unmasked genome Bowtie2 indexes:

```bash
mkdir -p fasta/index_bowtie_transcriptome fasta/index_bowtie_genome
bowtie2-build --threads 8 gencode.v47.pc_transcripts.fa fasta/index_bowtie_transcriptome/prot_coding_
bowtie2-build --threads 8 GRCh38.primary_assembly.genome.fa fasta/index_bowtie_genome/GRCh38_primary_genome
```

The original masked genome used an [AT-stretch BED file](https://osf.io/x2zeu/files/osfstorage). Download it separately, then build the masked genome and STAR index:

```bash
bedtools maskfasta -fi GRCh38.primary_assembly.genome.fa \
  -bed AT_stretches_merged.bed -fo GRCh38.primary_masked_AT.fa
mkdir -p fasta/index_masked_genome_star
bash scripts/star.sh --runMode genomeGenerate --runThreadN 8 \
  --genomeDir fasta/index_masked_genome_star \
  --genomeFastaFiles GRCh38.primary_masked_AT.fa \
  --sjdbGTFfile gencode.v47.primary_assembly.basic.annotation.gtf \
  --sjdbOverhang 79
```

Adjust `--sjdbOverhang` to read length minus one and thread counts to your machine. Generating and running the full human STAR index needs substantial RAM and disk space. Keep reference files outside the repository if convenient; absolute paths in `config.yaml` work.

## Analyse a sample

Edit `config.yaml` with a unique `sample` name, paired `R1` and `R2` FASTQ paths (plain `.fq`/`.fastq` or `.gz`), the Bowtie2 index **prefixes**, the STAR index **directory**, the annotation RDS path, and `threads`. Input paths can be relative to the repository root or absolute. R1 and R2 must have matching read names in the same order. The R2 sequence must begin with a 15-base UMI and `GTCAG`; reads shorter than 80 bases after UMI and adapter trimming are discarded.

```bash
snakemake --configfile config.yaml --cores 8
```

`--cores` is the machine-wide limit; `threads` is the per-alignment maximum. The three final files are `R_output/<sample>_final_big.rds`, `R_output/<sample>_final_small.rds`, and `R_output/<sample>_deduplication.csv`. The full table includes read-level alignment, tail composition, class, and type; the small table keeps the analysis columns without read identifiers, UMI, or raw sequence. The CSV reports counts before and after gene-plus-UMI deduplication.

Use a different sample name and FASTQ paths in a copied config file for another sample, for example `snakemake --configfile sample2.yaml --cores 8`. Intermediate files are sample-specific, so previous sample outputs can remain in place. For a dry run, add `-n`. If inputs or references change after a run, use `--forcerun` on the relevant rule or remove the affected outputs before rerunning.

## Troubleshooting

- **Missing reference index:** Check that `transcriptome_index` and `genome_index` are prefixes without `.1.bt2`/`.1.bt2l`, while `star_index` is a directory containing `Genome`.
- **No paired reads survived:** Inspect `fastq_read2/<sample>_fastp.html`. Confirm FASTQ pairing, R1 transcriptome matches, the R2 UMI/`GTCAG` motif, and read length.
- **R1 table conversion fails:** The transcriptome FASTA headers must include `transcript_id|gene_id`.
- **Annotation mismatch:** Gene IDs should be in the annotation RDS and chromosome names should match the STAR reference.
- **STAR reports zero input reads on Apple Silicon:** Run `bash scripts/install_star_macos.sh` in the active environment. The workflow checks for an empty SAM to catch this failure. See the [upstream issue](https://github.com/alexdobin/STAR/issues/2632).
- **STAR memory failure:** Run the full human reference on a machine with sufficient RAM, or use the tiny demo first to check the installation.
