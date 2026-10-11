"""Paired-end AT-tail analysis. Run from the repository root."""

from pathlib import Path
from scripts.sample_sheet import load_samples

configfile: "config.yaml"

SAMPLES = load_samples(config)
SAMPLE_IDS = list(SAMPLES)
TRANSCRIPTOME_INDEX = str(config["transcriptome_index"])
GENOME_INDEX = str(config["genome_index"])
STAR_INDEX = str(config["star_index"])
ANNOTATION = str(config["annotation"])


def bowtie_first_index(prefix):
    for suffix in (".1.bt2", ".1.bt2l"):
        if Path(prefix + suffix).exists():
            return prefix + suffix
    return prefix + ".1.bt2"


THREADS = int(config.get("threads", 4))
if THREADS < 1:
    raise ValueError("config 'threads' must be at least 1")


rule all:
    input:
        expand("R_output/{sample}_final_big.rds", sample=SAMPLE_IDS),
        expand("R_output/{sample}_final_small.rds", sample=SAMPLE_IDS),
        expand("R_output/{sample}_deduplication.csv", sample=SAMPLE_IDS)


rule align_read1:
    input:
        fastq=lambda wc: SAMPLES[wc.sample]["R1"],
        index=bowtie_first_index(TRANSCRIPTOME_INDEX)
    output:
        "sam/{sample}_R1_aligned.sam"
    threads: THREADS
    params:
        index=TRANSCRIPTOME_INDEX
    shell:
        "bowtie2 --threads {threads} --no-unal -x {params.index:q} "
        "-U {input.fastq:q} -S {output:q}"


rule select_read_pairs:
    input:
        sam="sam/{sample}_R1_aligned.sam",
        r1=lambda wc: SAMPLES[wc.sample]["R1"],
        r2=lambda wc: SAMPLES[wc.sample]["R2"]
    output:
        r1="fastq_read1/{sample}_R1_selected.fq",
        r2="fastq_read2/{sample}_R2_selected.fq"
    shell:
        "python R_scripts/select_read_pairs.py --sam {input.sam:q} "
        "--r1 {input.r1:q} --r2 {input.r2:q} "
        "--out-r1 {output.r1:q} --out-r2 {output.r2:q}"


rule fastp_filtering:
    input:
        r1="fastq_read1/{sample}_R1_selected.fq",
        r2="fastq_read2/{sample}_R2_selected.fq"
    output:
        r1="fastq_read1/{sample}_R1_fp.fq",
        r2="fastq_read2/{sample}_R2_fp.fq",
        html="fastq_read2/{sample}_fastp.html",
        json="fastq_read2/{sample}_fastp.json"
    threads: THREADS
    shell:
        "fastp --thread {threads} -i {input.r1:q} -I {input.r2:q} "
        "-o {output.r1:q} -O {output.r2:q} -Q --umi --umi_loc=read2 "
        "--umi_len=15 --trim_front2=5 -l 80 "
        "-h {output.html:q} -j {output.json:q}"


rule align_read2_genome:
    input:
        fastq="fastq_read2/{sample}_R2_fp.fq",
        index=bowtie_first_index(GENOME_INDEX)
    output:
        sam="sam/{sample}_R2_genome.sam",
        tags="R_input/{sample}_R2_genome_mapped.tsv"
    threads: THREADS
    params:
        index=GENOME_INDEX
    shell:
        "bowtie2 --threads {threads} -k 2 --no-unal -x {params.index:q} "
        "-U {input.fastq:q} -S {output.sam:q} && "
        "python R_scripts/extract_read_tags.py -i {output.sam:q} -o {output.tags:q}"


rule align_read2_masked_genome:
    input:
        fastq="fastq_read2/{sample}_R2_fp.fq",
        index=str(Path(STAR_INDEX) / "Genome")
    output:
        "R2_aligned_to_masked_genome/{sample}/Aligned.out.sam"
    threads: THREADS
    params:
        index=STAR_INDEX,
        prefix=lambda wc: f"R2_aligned_to_masked_genome/{wc.sample}/"
    shell:
        "mkdir -p {params.prefix:q} && "
        "bash scripts/star.sh --genomeDir {params.index:q} --readFilesIn {input.fastq:q} "
        "--runThreadN {threads} --outFileNamePrefix {params.prefix:q} "
        "--outSAMtype SAM --outFilterMismatchNoverReadLmax 0.05 "
        "--alignIntronMax 50000 --outFilterMultimapNmax 1 "
        "--outSAMunmapped Within && "
        "python R_scripts/check_star_alignment.py {input.fastq:q} {output:q}"


rule extract_alignment_tables:
    input:
        r1="sam/{sample}_R1_aligned.sam",
        r2="R2_aligned_to_masked_genome/{sample}/Aligned.out.sam"
    output:
        r1="R_input/{sample}_R1_aligned.tsv",
        r2="R_input/{sample}_R2_aligned_to_masked_genome.tsv"
    shell:
        "python R_scripts/sam_to_tables.py --r1 {input.r1:q} --r2 {input.r2:q} "
        "--out-r1 {output.r1:q} --out-r2 {output.r2:q}"


rule tail_analysis:
    input:
        tabR1="R_input/{sample}_R1_aligned.tsv",
        tabR2="R_input/{sample}_R2_aligned_to_masked_genome.tsv",
        tabR2_2="R_input/{sample}_R2_genome_mapped.tsv",
        annotation=ANNOTATION
    output:
        big_table="R_output/{sample}_final_big.rds",
        small_table="R_output/{sample}_final_small.rds",
        deduplication="R_output/{sample}_deduplication.csv"
    script:
        "R_scripts/run_tail_analysis.R"
