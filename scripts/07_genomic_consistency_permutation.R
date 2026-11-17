# Reproduce the manuscript's CASE-LEVEL cluster-label permutation test.
# NEVER export sequence names, full sequences, individual case labels or pairs.
# Both files MUST be prepared inside the authorized private environment:
#   genomic_cluster_labels_private.csv: genomic_index, cluster_id (1 row/case)
#   genomic_pairwise_snp_private.csv: index_A, index_B, SNP_distance
# Genomic index is 1..N in exactly the same order as indices in the pairs file;
# include only sequenced cases belonging to multi-case reconstructed clusters.
# Original Rmd default minimum_pairwise_comparable_fraction=NULL (no pair QC).
source("R/common.R")
source("R/genomics_core.R")
labels <- read_private_csv("genomic_cluster_labels_private.csv")
require_columns(labels,c("genomic_index","cluster_id"),"genomic cluster labels")
if (anyNA(labels$cluster_id) || anyNA(labels$genomic_index) ||
    any(labels$genomic_index != as.integer(labels$genomic_index)) ||
    anyDuplicated(labels$genomic_index) ||
    !identical(as.integer(labels$genomic_index),seq_len(nrow(labels)))) {
  stop("genomic_index must be ordered 1..number of clustered sequenced cases")
}
pair <- read_private_csv("genomic_pairwise_snp_private.csv")
res <- case_label_genomic_permutation(labels$cluster_id,pair,
  thresholds=c(1L,2L,3L,5L),permutations=9999L,seed=20260807L)
safe_write_csv(res,"genomic_case_level_permutation_summary_PRIVATE.csv")
message("Manuscript two-sided SNP permutation test completed; only private aggregate results written.")
