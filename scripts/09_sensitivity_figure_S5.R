# Figure S5 plotting data from the final Rmd's sensitivity metrics.
# These are descriptive relative changes vs the case-major baseline scenario.
source("R/common.R")
if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2")
f <- private_result_path("network_sensitivity_results_numeric.csv")
if (!file.exists(f)) stop("Run scripts/03_parameter_sensitivity.R first")
d <- readr::read_csv(f,show_col_types=FALSE)
metrics <- c("mean_Ri","links_within_500m_pct","unconnected_cases_pct","duration_diameter_slope")
require_columns(d,c("parameter_group","parameter_setting","parameter_order","is_baseline",metrics),"sensitivity")
b <- d[which(d$is_baseline),,drop=FALSE]
if (nrow(b)!=1L) stop("One baseline is required")
groups <- unique(d$parameter_group[!d$is_baseline])
plot_data <- do.call(rbind,lapply(groups,function(g) {
  sub <- d[d$parameter_group==g,,drop=FALSE]
  baseline_row <- b
  baseline_row$parameter_group <- g
  baseline_row$parameter_setting <- switch(g,
    "Serial interval mean"="8","Serial interval SD"="3",
    "Spatial decay coefficient"="1/250","Maximum candidate distance"="1300",
    "Background source weight"="10^-4")
  baseline_row$parameter_order <- switch(g,
    "Serial interval mean"=3L,"Serial interval SD"=2L,
    "Spatial decay coefficient"=5L,"Maximum candidate distance"=2L,
    "Background source weight"=2L)
  sub <- rbind(sub,baseline_row)
  sub <- sub[order(sub$parameter_order),,drop=FALSE]
  do.call(rbind,lapply(metrics,function(metric) {
    v0 <- b[[metric]][1]
    if (!is.finite(v0) || abs(v0)<1e-10) stop("Invalid baseline for ",metric)
    data.frame(parameter_group=g,parameter_setting=sub$parameter_setting,
      parameter_order=sub$parameter_order,metric=metric,
      relative_change_pct=100*(sub[[metric]]-v0)/abs(v0))
  }))
}))
plot_data$parameter_group <- factor(plot_data$parameter_group,levels=groups)
plt <- ggplot2::ggplot(plot_data,ggplot2::aes(x=parameter_order,y=relative_change_pct,
    group=metric,color=metric)) +
  ggplot2::geom_hline(yintercept=0,linetype="dashed",color="grey65") +
  ggplot2::geom_line(linewidth=.7) + ggplot2::geom_point(size=1.5) +
  ggplot2::facet_wrap(~parameter_group,scales="free_x",ncol=2) +
  ggplot2::labs(x="Parameter setting (ordered; see plotted data)",
    y="Relative change from scenario baseline (%)",color="Metric") +
  ggplot2::theme_bw()
ggplot2::ggsave(private_result_path("Figure_S5_sensitivity_profiles.png"),plt,width=12,height=8,dpi=300)
safe_write_csv(plot_data,"Figure_S5_aggregated_plotting_data.csv")
