# Offline base-R synthetic invariants. From project root: Rscript tests/test_core.R
source("R/common.R")
source("R/network_core.R")
source("R/cluster_core.R")
source("R/knox_core.R")
source("R/genomics_core.R")
small <- data.frame(case_id=c("TEST_A","TEST_B","TEST_C","TEST_D"),
                    onset_date=c("2025-07-01","2025-07-09","2025-07-11","2025-07-15"),
                    case_type=c("Imported","Local","Local","Imported"),
                    utm_x=c(500000,500100,500900,500200),
                    utm_y=c(2600000,2600100,2600900,2600300))
d <- validate_network_input(small)
p <- network_parameters(); p$n_samples <- 300L
scores <- make_score_matrix(d,p)
stopifnot(all(scores[1,]==0),all(scores[4,]==0),scores[2,1]>0,
          scores[2,3]==0,scores[3,2]>0)
draws <- sample_networks(d,scores,p,rng_order="simulation_major")
stopifnot(identical(dim(draws),c(4L,300L)),all(draws[c(1,4),]==0L))
stopifnot(all(draws[2,] %in% c(0L,1L)),all(draws[3,] %in% c(0L,1L,2L)))
# Determinism of original simulation-major draw order.
stopifnot(identical(draws,sample_networks(d,scores,p,rng_order="simulation_major")))
case_major <- sample_networks(d,scores,p,rng_order="case_major")
stopifnot(identical(case_major,sample_networks(d,scores,p,rng_order="case_major")))
res <- extract_network_summary(d,draws)
stopifnot(abs(sum(res$mean_Ri)-sum(draws>0)/ncol(draws))<1e-10)
# The background state must win when most frequent (even if observed candidate exists).
# Also check tie: 0 wins; other ties resolve to smallest NUMERIC index.
trial <- matrix(0L,nrow=12,ncol=6)
trial[1,] <- c(0,0,0,0,2,2)
trial[2,] <- c(2,2,10,10,0,0)
modal <- posterior_source_summary(trial)
stopifnot(modal$parent[1]==0L, modal$parent[2]==0L,
          modal$n_tied_sources[2]==3L, modal$support_prob[1]==4/6)
roots <- find_cluster_root(c(0L,1L,2L,0L))
stopifnot(identical(roots,c(1L,1L,1L,4L)))
# A seed with no upstream but one downstream must NOT be called unconnected.
network <- d
network$top_parent_id <- c(0L,1L,2L,0L)
cl <- cluster_statistics(network)
stopifnot(nrow(cl$clusters)==1L,cl$clusters$cluster_size==3L,
          cl$clusters$maximum_generation_depth==2L,
          sum(cl$membership$is_unconnected)==1L,
          sum(cl$membership$generation_depth==0L)==2L)
knox <- knox_grid(as.matrix(d[,c("utm_x","utm_y")]),d$onset_date,
                  spatial_m=c(100,1500),temporal_days=7,permutations=19L,seed=123L)
stopifnot(nrow(knox)==2L,all(knox$permutation_p>0),all(knox$permutation_p<=1))
# Genomic tests: one row per pair; several dependent pairs share a case.
idx <- t(utils::combn(5,2))
pairwise <- data.frame(index_A=idx[,1],index_B=idx[,2],
                       SNP_distance=c(0L,1L,2L,3L,1L,1L,2L,4L,5L,0L))
g <- case_label_genomic_permutation(c(1,1,2,2,2),pairwise,permutations=39L,seed=101L)
stopifnot(nrow(g)==4L,identical(g$SNP_threshold,c(1L,2L,3L,5L)),
          all(g$permutation_p_value >=1/40),all(g$permutation_p_value<=1),
          all(g$holm_adjusted_permutation_p>=g$permutation_p_value),
          identical(g,case_label_genomic_permutation(c(1,1,2,2,2),pairwise,
                                                     permutations=39L,seed=101L)))
cat("PASS: synthetic reconstruction, background mode, forest clusters, Knox, genomic permutation\n")
