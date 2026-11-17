# Analytic figure COMPONENTS (not certified manuscript layouts).
# Use only approved restricted derived data and keep all graphics git-ignored.
source("R/common.R")
source("R/cluster_core.R")
if (!requireNamespace("ggplot2",quietly=TRUE) || !requireNamespace("scales",quietly=TRUE)) {
  stop("Install ggplot2 and scales")
}
plot_dir <- dirname(private_result_path("figure_placeholder.png"))
save_plot <- function(p,name,w=7,h=5) {
  ggplot2::ggsave(file.path(plot_dir,name),p,width=w,height=h,dpi=300,bg="white")
}
network_path <- private_result_path("network_case_summary.csv")
if (file.exists(network_path)) {
  d <- readr::read_csv(network_path,show_col_types=FALSE)
  require_columns(d,c("id_int","case_id","onset_date","mean_Ri","Ri_LCI","Ri_UCI",
                      "top_parent_id","infect_distance_m"),"network case summary")
  d$onset_date <- as.Date(d$onset_date)
  # Figure 2A: mirror the original Rmd's grouped-on-onset-date 7-row right
  # moving mean; "Rt" in the legacy plot was NOT incidence-based R_t.
  days <- split(seq_len(nrow(d)), d$onset_date)
  daily <- do.call(rbind,lapply(days,function(idx) data.frame(
    onset_date=d$onset_date[idx[1]], model_mean=mean(d$mean_Ri[idx]),
    lower_individual_mean=mean(d$Ri_LCI[idx]),
    upper_individual_mean=mean(d$Ri_UCI[idx]), case_count=length(idx))))
  daily <- daily[order(daily$onset_date),,drop=FALSE]
  rolling7 <- function(x) as.numeric(stats::filter(x,rep(1/7,7),sides=1))
  daily$smoothed <- rolling7(daily$model_mean)
  daily$lower_smooth <- rolling7(daily$lower_individual_mean)
  daily$upper_smooth <- rolling7(daily$upper_individual_mean)
  g2a <- ggplot2::ggplot(daily,ggplot2::aes(x=onset_date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower_smooth,ymax=upper_smooth),alpha=.15) +
    ggplot2::geom_line(ggplot2::aes(y=smoothed),linewidth=.8) +
    ggplot2::labs(x="Symptom-onset date",y="7-observed-date smoothed mean inferred direct offspring") +
    ggplot2::theme_classic()
  save_plot(g2a,"Figure_2A_component.png")
  # Figure 2B: cumulative observed-source residential link distances.
  dist <- d$infect_distance_m[is.finite(d$infect_distance_m)]
  if (length(dist)) {
    links <- data.frame(distance_m=dist)
    g2b <- ggplot2::ggplot(links,ggplot2::aes(distance_m)) +
      ggplot2::stat_ecdf(geom="step",linewidth=.85) +
      ggplot2::geom_vline(xintercept=500,linetype="dashed") +
      ggplot2::labs(x="Residential link distance (m)",y="Cumulative fraction of inferred links") +
      ggplot2::theme_classic()
    save_plot(g2b,"Figure_2B_component.png")
  }
  # Figure 2C: same >=2-case clusters and unfiltered duration-diameter OLS
  # used in the reported descriptive slope.
  cl <- cluster_statistics(d)$clusters
  if (nrow(cl)>=3L) {
    g2c <- ggplot2::ggplot(cl,ggplot2::aes(duration_days,residential_diameter_m)) +
      ggplot2::geom_smooth(method="lm",formula=y~x,se=TRUE,colour="grey25") +
      ggplot2::geom_point(ggplot2::aes(size=cluster_size),alpha=.7) +
      ggplot2::labs(x="Cluster duration (days)",y="Residential cluster diameter (m)",
                    size="Cluster size") + ggplot2::theme_classic()
    save_plot(g2c,"Figure_2C_component.png")
  }
  # Figure 2D: the original top-20% share is derived from the full parent
  # matrix (all 10,000 networks), NOT the modal descriptive network.
  parents_file <- private_result_path("sim_parents_10000.rds")
  if (file.exists(parents_file)) {
    parents <- readRDS(parents_file)
    offspring <- tabulate(as.integer(parents[parents>0L]),nbins=nrow(parents))/ncol(parents)
    decreasing <- sort(offspring,decreasing=TRUE)
    share <- sum(decreasing[seq_len(ceiling(.2*length(decreasing)))])/sum(decreasing)
    curve <- data.frame(proportion_cases=c(0,seq_along(decreasing)/length(decreasing)),
                        proportion_onward=c(0,cumsum(decreasing)/sum(decreasing)))
    if (!all(is.finite(curve$proportion_onward))) stop("No inferred offspring")
    g2d <- ggplot2::ggplot(curve,ggplot2::aes(proportion_cases,proportion_onward)) +
      ggplot2::geom_area(fill="#DCE6F2",alpha=.45) +
      ggplot2::geom_line(color="#3E648D",linewidth=1.1) +
      ggplot2::annotate("segment",x=.2,xend=.2,y=0,yend=share,linetype="dashed",color="#C55A5A") +
      ggplot2::scale_x_continuous(labels=scales::label_percent()) +
      ggplot2::scale_y_continuous(labels=scales::label_percent()) +
      ggplot2::labs(x="Top fraction of cases ranked by inferred direct offspring",
                    y="Cumulative share of inferred onward transmission") + ggplot2::theme_classic()
    save_plot(g2d,"Figure_2D_component.png")
  }
}
# Figure 3A (Model 3 linear/categorical coefficients) and 3B/C (marginal
# standardized LATENT predictor). These are source-data components only.
coef_path <- private_result_path("primary_coefficients.csv")
if (file.exists(coef_path)) {
  cf <- readr::read_csv(coef_path,show_col_types=FALSE)
  cf <- subset(cf,model=="M3" & term %in% c("case_typeLocal","z_BI","z_outflow",
                                            "z_dens_park","z_log_hosp","z_log_trans"))
  if (nrow(cf)) {
    g3a <- ggplot2::ggplot(cf,ggplot2::aes(x=estimate,y=stats::reorder(term,estimate))) +
      ggplot2::geom_vline(xintercept=0,linetype="dashed") +
      ggplot2::geom_segment(ggplot2::aes(x=ci_lower,xend=ci_upper,
        y=stats::reorder(term,estimate),yend=stats::reorder(term,estimate))) +
      ggplot2::geom_point(size=2) +
      ggplot2::labs(x="Model 3 Tobit coefficient (95% Wald CI)",y=NULL) +
      ggplot2::theme_classic()
    save_plot(g3a,"Figure_3A_component.png")
  }
}
for (target in c("onset_day","log_dist_market")) {
  f <- private_result_path(paste0("figure3_",target,"_latent_curve.csv"))
  if (!file.exists(f)) next
  d <- readr::read_csv(f,show_col_types=FALSE)
  g <- ggplot2::ggplot(d,ggplot2::aes(x,latent_mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=.18) +
    ggplot2::geom_line(linewidth=.85) +
    ggplot2::labs(x=if (target=="onset_day") "Days since first onset" else "Log distance to market",
                  y="Adjusted latent Tobit prediction") + ggplot2::theme_classic()
  save_plot(g,paste0("Figure_3_",target,"_component.png"))
}
# Figure S4: distribution of linear coefficients and joint-spline significance.
cf <- private_result_path("sampled_network_all_coefficients.csv")
sp <- private_result_path("sampled_network_spline_stability.csv")
if (file.exists(cf)) {
  d <- readr::read_csv(cf,show_col_types=FALSE)
  d <- subset(d,term %in% c("case_typeLocal","z_BI","z_outflow",
                           "z_dens_park","z_log_hosp","z_log_trans"))
  g <- ggplot2::ggplot(d,ggplot2::aes(estimate,term)) +
    ggplot2::geom_vline(xintercept=0,linetype="dashed",colour="grey50") +
    ggplot2::geom_violin(fill="#3B6FA8",alpha=.25) +
    ggplot2::geom_boxplot(width=.2,outlier.shape=NA,alpha=.55) +
    ggplot2::labs(x="Tobit coefficient across sampled networks",y=NULL) +
    ggplot2::theme_classic()
  save_plot(g,"Figure_S4A_component.png")
}
if (file.exists(sp)) {
  d <- readr::read_csv(sp,show_col_types=FALSE)
  g <- ggplot2::ggplot(d,ggplot2::aes(p_less_0_05_percent,
      stats::reorder(variable,p_less_0_05_percent))) +
    ggplot2::geom_col(fill="#79589F") +
    ggplot2::scale_x_continuous(limits=c(0,100)) +
    ggplot2::labs(x="Percentage of models with spline joint Wald P < 0.05",y=NULL) +
    ggplot2::theme_classic()
  save_plot(g,"Figure_S4B_component.png")
}
message("Available analytic figure components generated under results/private.")
