# Final Rmd's one-at-a-time network sensitivity, not the exploratory 1,000-draw
# version. IMPORTANT: 10,000 samples PER scenario using the late Rmd's distinct
# case-major draw order. Includes source 0 in modal assignment and isolates.
source("R/common.R")
source("R/network_core.R")
source("R/cluster_core.R")
cases <- validate_network_input(read_private_csv("cases_network.csv"))
base <- network_parameters()
scenario <- list(list(group="Baseline", setting="Baseline", order=1L,
                      settings=list(), baseline=TRUE))
for (i in seq_along(c(6,7,9,10))) {
  x <- c(6,7,9,10)[i]
  scenario[[length(scenario)+1L]] <- list(group="Serial interval mean",setting=as.character(x),
    order=c(1L,2L,4L,5L)[i],settings=list(si_mean=x),baseline=FALSE)
}
for (i in seq_along(c(2,4,6))) {
  x <- c(2,4,6)[i]
  scenario[[length(scenario)+1L]] <- list(group="Serial interval SD",setting=as.character(x),
    order=c(1L,3L,4L)[i],settings=list(si_sd=x),baseline=FALSE)
}
for (i in seq_along(c(50,100,150,200,300,400,500,1000))) {
  x <- c(50,100,150,200,300,400,500,1000)[i]
  scenario[[length(scenario)+1L]] <- list(group="Spatial decay coefficient",
    setting=paste0("1/",x),order=c(1L,2L,3L,4L,6L,7L,8L,9L)[i],
    settings=list(distance_scale_m=x),baseline=FALSE)
}
for (i in seq_along(c(500,2000,3000))) {
  x <- c(500,2000,3000)[i]
  scenario[[length(scenario)+1L]] <- list(group="Maximum candidate distance",
    setting=as.character(x),order=c(1L,3L,4L)[i],
    settings=list(distance_cutoff_m=x),baseline=FALSE)
}
for (i in seq_along(c(1e-5,1e-3))) {
  x <- c(1e-5,1e-3)[i]
  scenario[[length(scenario)+1L]] <- list(group="Background source weight",
    setting=c("10^-5","10^-3")[i],order=c(1L,3L)[i],
    settings=list(background_weight=x),baseline=FALSE)
}
run_scenario <- function(sc) {
  par <- modifyList(base,sc$settings)
  if (par$n_samples != 10000L) stop("Paper sensitivity uses exactly 10,000 draws")
  message("Sensitivity: ",sc$group," = ",sc$setting)
  draws <- sample_networks(cases,make_score_matrix(cases,par),par,rng_order="case_major")
  net <- extract_network_summary(cases,draws,max_distance_m=par$distance_cutoff_m)
  cl <- cluster_statistics(net)
  slope <- cluster_slope(cl$clusters)
  distances <- net$infect_distance_m[is.finite(net$infect_distance_m)]
  onset_lags <- net$infect_time_diff_d[is.finite(net$infect_time_diff_d)]
  n_links <- length(distances)
  n <- nrow(net)
  n_unconnected <- sum(net$is_unconnected)
  ci_500 <- wilson_ci(sum(distances<=500),n_links)
  ci_1000 <- wilson_ci(sum(distances<=1000),n_links)
  ci_unconnected <- wilson_ci(n_unconnected,n)
  # Source calculates a per-network mean direct-offspring count and takes its
  # 2.5th and 97.5th percentiles. This is not a case-level CI.
  mean_ri_per_draw <- colSums(draws>0L)/n
  mean_ri_ci <- as.numeric(stats::quantile(mean_ri_per_draw,c(.025,.975),names=FALSE))
  lag_ci <- if (length(onset_lags)) as.numeric(stats::quantile(onset_lags,c(.025,.975),names=FALSE)) else c(NA_real_,NA_real_)
  out <- data.frame(scenario_name=if (sc$baseline) "Baseline" else paste(sc$group,sc$setting,sep="_"),
    parameter_group=sc$group,parameter_setting=sc$setting,parameter_order=sc$order,
    is_baseline=sc$baseline,
    theta=1/par$distance_scale_m,maximum_candidate_distance_m=par$distance_cutoff_m,
    background_source_weight=par$background_weight,serial_interval_mean_days=par$si_mean,
    serial_interval_sd_days=par$si_sd,n_sampled_networks=ncol(draws),
    mean_Ri=mean(mean_ri_per_draw),mean_Ri_lower=mean_ri_ci[1],mean_Ri_upper=mean_ri_ci[2],
    n_observed_source_links=n_links,
    median_link_distance_m=if (n_links) median(distances) else NA_real_,
    link_distance_q1_m=if (n_links) unname(stats::quantile(distances,.25)) else NA_real_,
    link_distance_q3_m=if (n_links) unname(stats::quantile(distances,.75)) else NA_real_,
    links_within_500m_pct=if (n_links) 100*mean(distances<=500) else NA_real_,
    links_within_500m_lower=100*ci_500["lower"],links_within_500m_upper=100*ci_500["upper"],
    links_within_1000m_pct=if (n_links) 100*mean(distances<=1000) else NA_real_,
    links_within_1000m_lower=100*ci_1000["lower"],links_within_1000m_upper=100*ci_1000["upper"],
    n_clusters=nrow(cl$clusters),
    duration_diameter_slope=slope$estimate,
    duration_diameter_slope_lower=slope$lower,duration_diameter_slope_upper=slope$upper,
    duration_diameter_slope_p=slope$p_value,duration_diameter_r2=slope$r_squared,
    unconnected_cases_n=n_unconnected,unconnected_cases_pct=100*n_unconnected/n,
    unconnected_cases_lower=100*ci_unconnected["lower"],
    unconnected_cases_upper=100*ci_unconnected["upper"],
    cases_without_upstream_n=sum(!net$has_observed_upstream),
    locally_acquired_background_n=sum(net$case_type=="local" & net$top_parent_id==0L),
    locally_acquired_background_pct=100*sum(net$case_type=="local" & net$top_parent_id==0L)/sum(net$case_type=="local"),
    selected_link_onset_interval_mean=if (length(onset_lags)) mean(onset_lags) else NA_real_,
    selected_link_onset_interval_lower=lag_ci[1],selected_link_onset_interval_upper=lag_ci[2],
    maximum_source_tie_cases=sum(net$n_tied_sources>1L))
  out
}
result <- do.call(rbind,lapply(scenario,run_scenario))
safe_write_csv(result,"network_sensitivity_results_numeric.csv")
message("21 sensitivity scenarios completed; outputs restricted to results/private.")
