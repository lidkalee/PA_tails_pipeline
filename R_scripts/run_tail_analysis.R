library(stringr)
library(stringi)
library(dplyr)
source("R_scripts/run_AT_tail_analysis_df.R")
source("R_scripts/assign_type_and_class.R")

read_tsv <- function(path, columns) {
  if (file.info(path)$size == 0) {
    names <- paste0("V", seq_len(columns))
    return(stats::setNames(as.data.frame(matrix(character(), nrow = 0, ncol = columns)), names))
  }
  read.table(path, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
}

R1 <- read_tsv(snakemake@input[["tabR1"]], 7)
R2 <- read_tsv(snakemake@input[["tabR2"]], 6)
if (nrow(R1) == 0) stop("No R1 alignments: check the FASTQ and transcriptome index")

R2$UMI <- stringr::str_extract(R2$V1, "(?<=:)[ACGTN]{15}$")
R2$V1 <- sub(":(?:[ACGTN]{15})$", "", R2$V1)
merged_dat <- dplyr::inner_join(
  dplyr::select(R1, V1, V2, V3, V4, V5, V7),
  dplyr::select(R2, V1, V2, V3, V4, V5, V6, UMI),
  by = "V1"
)
if (nrow(merged_dat) == 0) stop("No paired reads survived alignment and filtering; inspect fastp reports")
colnames(merged_dat) <- c("read_name", "geneR1", "transR1", "flagR1", "mapR1", "CIGAR_R1",
                          "flagR2", "chrR2", "posR2", "CIGAR_R2", "seq_R2", "UMI")
if (anyNA(merged_dat$UMI)) stop("UMI missing from fastp read names")

dedup_dat <- dplyr::distinct(merged_dat, geneR1, UMI, .keep_all = TRUE)
dedup_dat$geneR1 <- sub("\\..*", "", dedup_dat$geneR1)
dedup_stats <- data.frame(stage = c("before_dedup", "after_dedup"),
                          reads = c(nrow(merged_dat), nrow(dedup_dat)))
annotation <- readRDS(snakemake@input[["annotation"]])
required <- c("gene_id", "seqnames", "start", "end", "strand")
if (!all(required %in% names(annotation))) stop("Annotation RDS lacks required columns")
dedup_dat <- dplyr::left_join(dedup_dat, annotation[, required], by = c("geneR1" = "gene_id"))
dedup_dat <- run_AT_tail_analysis_df(dedup_dat)
dedup_dat <- assign_type_and_class(dedup_dat)

R2_genome <- read_tsv(snakemake@input[["tabR2_2"]], 2)
R2_genome$V1 <- sub(":(?:[ACGTN]{15})$", "", R2_genome$V1)
unique_genome <- R2_genome$V1[R2_genome$V2 == "unique"]
dedup_dat$al_ete_bt <- ifelse(dedup_dat$tail_len < 40 & dedup_dat$type == "grep" &
                              dedup_dat$read_name %in% unique_genome, "yes", "no")

small_columns <- c("geneR1", "flagR1", "flagR2", "chrR2", "posR2", "CIGAR_R2",
                   "seqnames", "start", "end", "strand", "A_count", "T_count",
                   "G_count", "C_len", "tail_len", "high_AT_mix", "tail",
                   "type", "class", "al_ete_bt")
saveRDS(dedup_dat, snakemake@output[["big_table"]])
saveRDS(dedup_dat[, small_columns, drop = FALSE], snakemake@output[["small_table"]])
write.csv(dedup_stats, snakemake@output[["deduplication"]], row.names = FALSE)
