# Case-level cluster-label permutation for SNP similarity (late Rmd, section 23.1).
# Input is an OFFLINE RESTRICTED analysis representation, not FASTA/name linkage:
# a cluster label per sequenced case and pairs of integer indices plus SNP distance.
# This is the manuscript's two-sided empirical permutation test, not a chi-square
# test on pair counts (pairs share cases and are not independent observations).
case_label_genomic_permutation <- function(cluster_labels, pairwise,
                                           thresholds=c(1L,2L,3L,5L),
                                           permutations=9999L,
                                           seed=20260807L) {
  n <- length(cluster_labels)
  if (n < 3L || anyNA(cluster_labels)) stop("Require >=3 nonmissing case-level cluster labels")
  if (!is.data.frame(pairwise) || !all(c("index_A","index_B","SNP_distance") %in% names(pairwise))) {
    stop("pairwise must have index_A, index_B, SNP_distance")
  }
  ia <- pairwise$index_A
  ib <- pairwise$index_B
  distance <- pairwise$SNP_distance
  if (!length(ia) || anyNA(ia) || anyNA(ib) || anyNA(distance) ||
      any(ia != as.integer(ia)) || any(ib != as.integer(ib)) ||
      any(ia < 1 | ia > n | ib < 1 | ib > n | ia >= ib) ||
      any(!is.finite(distance)) || any(distance < 0) ||
      any(distance != as.integer(distance))) stop("Invalid pair indices or SNP distances")
  # The analysis uses each undirected pair at most once.
  pair_keys <- paste(ia,ib,sep="-")
  if (anyDuplicated(pair_keys)) stop("Repeated pairwise comparisons")
  if (length(thresholds) != 4L || !identical(as.numeric(thresholds),c(1,2,3,5))) {
    stop("Published SNP thresholds must be c(1,2,3,5)")
  }
  if (!is.numeric(permutations) || length(permutations)!=1L ||
      is.na(permutations) || permutations<1L || permutations!=as.integer(permutations)) {
    stop("permutations must be a positive integer")
  }
  threshold_indicator <- vapply(thresholds,
    function(t) distance <= t, logical(length(distance)))
  threshold_indicator <- matrix(threshold_indicator, nrow=length(distance), ncol=length(thresholds))
  statistic <- function(labels, with_counts=FALSE) {
    within <- labels[ia] == labels[ib]
    n_within <- sum(within)
    n_between <- length(within)-n_within
    if (n_within==0L || n_between==0L) stop("No valid within- and between-cluster comparison")
    similar_within <- colSums(threshold_indicator[within,,drop=FALSE])
    similar_between <- colSums(threshold_indicator[!within,,drop=FALSE])
    diff <- similar_within/n_within-similar_between/n_between
    if (with_counts) {
      list(diff=diff, similar_within=similar_within,
           similar_between=similar_between,n_within=n_within,n_between=n_between)
    } else diff
  }
  obs <- statistic(cluster_labels,with_counts=TRUE)
  set.seed(as.integer(seed))
  perm_stats <- matrix(NA_real_,nrow=as.integer(permutations),ncol=length(thresholds))
  for (i in seq_len(as.integer(permutations))) {
    perm_stats[i,] <- statistic(sample(cluster_labels,replace=FALSE))
  }
  two_p <- (colSums(abs(perm_stats) >= matrix(abs(obs$diff), nrow=nrow(perm_stats),
                               ncol=length(thresholds),byrow=TRUE))+1)/(permutations+1)
  one_p <- (colSums(perm_stats >= matrix(obs$diff,nrow=nrow(perm_stats),
                               ncol=length(thresholds),byrow=TRUE))+1)/(permutations+1)
  lower <- apply(perm_stats,2,stats::quantile,probs=.025,names=FALSE)
  median <- apply(perm_stats,2,stats::median)
  upper <- apply(perm_stats,2,stats::quantile,probs=.975,names=FALSE)
  data.frame(SNP_threshold=thresholds,
             within_similar_n=as.integer(obs$similar_within),
             within_total_n=obs$n_within,
             within_similar_percent=100*obs$similar_within/obs$n_within,
             between_similar_n=as.integer(obs$similar_between),
             between_total_n=obs$n_between,
             between_similar_percent=100*obs$similar_between/obs$n_between,
             difference_percentage_points=100*obs$diff,
             null_difference_2.5th_percentile=100*lower,
             null_difference_median=100*median,
             null_difference_97.5th_percentile=100*upper,
             permutation_p_value=two_p,
             permutation_p_value_one_sided=one_p,
             holm_adjusted_permutation_p=stats::p.adjust(two_p,method="holm"))
}
