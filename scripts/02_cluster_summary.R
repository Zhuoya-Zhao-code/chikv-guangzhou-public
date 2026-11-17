source("R/common.R")
source("R/cluster_core.R")
d <- readr::read_csv(private_result_path("network_case_summary.csv"),show_col_types=FALSE)
network <- cluster_statistics(d)
safe_write_csv(network$membership,"cluster_membership_PRIVATE.csv")
safe_write_csv(network$clusters,"cluster_characteristics_PRIVATE.csv")
slope <- cluster_slope(network$clusters)
summary <- data.frame(n_clusters=nrow(network$clusters),
                      n_connected_cases=sum(network$clusters$cluster_size),
                      n_unconnected=sum(network$membership$is_unconnected),
                      largest_cluster=max(network$clusters$cluster_size),
                      maximum_generation_depth=max(network$membership$generation_depth),
                      diameter_duration_slope_m_day=slope$estimate,
                      diameter_duration_slope_lower=slope$lower,
                      diameter_duration_slope_upper=slope$upper,
                      diameter_duration_slope_p=slope$p_value,
                      diameter_duration_r_squared=slope$r_squared)
safe_write_csv(summary,"cluster_aggregate_summary.csv")
message("Cluster metrics use all multi-case components, including duration=0 when present.")
