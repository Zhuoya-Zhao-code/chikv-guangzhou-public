# Primary reconstruction: preserve simulation-major draws and late modal
# source selection INCLUDING background 0. Outputs are strictly private.
source("R/common.R")
source("R/network_core.R")
source("R/cluster_core.R")
params <- network_parameters()
cases <- validate_network_input(read_private_csv("cases_network.csv"))
parents <- sample_networks(cases, make_score_matrix(cases, params), params,
                           rng_order="simulation_major")
summary_private <- extract_network_summary(cases, parents, max_distance_m=params$distance_cutoff_m)
clusters <- cluster_statistics(summary_private)
safe_write_csv(summary_private, "network_case_summary.csv")
safe_write_csv(summary_private[, c("case_id", "id_int")], "network_case_map.csv")
saveRDS(parents, private_result_path("sim_parents_10000.rds"), compress=FALSE)
link_dist <- summary_private$infect_distance_m[is.finite(summary_private$infect_distance_m)]
if (!length(link_dist)) stop("No observed-source links selected")
ci500 <- wilson_ci(sum(link_dist<=500),length(link_dist))
ci_unconnected <- wilson_ci(sum(summary_private$is_unconnected),nrow(summary_private))
ri_decreasing <- sort(summary_private$mean_Ri,decreasing=TRUE)
res <- data.frame(
  n_cases=nrow(cases), n_sampled_networks=ncol(parents),
  n_imported=sum(cases$case_type=="imported"),
  n_observed_source_links=sum(summary_private$top_parent_id!=0L),
  n_clusters=nrow(clusters$clusters), n_connected_cases=sum(clusters$clusters$cluster_size),
  n_unconnected=sum(summary_private$is_unconnected),
  pct_unconnected=100*mean(summary_private$is_unconnected),
  unconnected_wilson_lower_pct=100*ci_unconnected["lower"],
  unconnected_wilson_upper_pct=100*ci_unconnected["upper"],
  largest_cluster=max(clusters$clusters$cluster_size),
  max_generation_depth=max(clusters$membership$generation_depth),
  median_link_distance_m=median(link_dist),
  link_500_pct=100*mean(link_dist<=500),
  link_500_wilson_lower_pct=100*ci500["lower"],
  link_500_wilson_upper_pct=100*ci500["upper"],
  top20pct_onward_share=sum(ri_decreasing[seq_len(ceiling(.2*length(ri_decreasing)))])/sum(ri_decreasing),
  n_modal_tied_cases=sum(summary_private$n_tied_sources>1L)
)
safe_write_csv(res,"network_aggregate_summary.csv")
message("Primary reconstruction complete; all case-level artifacts remain git-ignored.")
