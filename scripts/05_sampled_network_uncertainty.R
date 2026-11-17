# Final 300 sampled-network Tobit uncertainty (case-level mappings remain private).
# Based on chikv_final_remaining_analysis(1).R, NOT the placeholder CSV variant.
source("R/common.R")
source("R/regression_core.R")
if (!requireNamespace("AER", quietly = TRUE)) stop("Install R package AER")
d <- prepare_regression(read_private_csv("Case_Data_Phase1_Ready.csv"))
require_columns(d, "case_id", "Case_Data_Phase1_Ready.csv")
check_case_ids(d$case_id)
parent_path <- private_result_path("sim_parents_10000.rds")
if (!file.exists(parent_path)) stop("First run scripts/01_reconstruct.R")
parents <- readRDS(parent_path)
if (!is.matrix(parents) || ncol(parents) < 300L || anyNA(parents) ||
    any(parents < 0L | parents > nrow(parents)) || any(parents != as.integer(parents))) {
  stop("Invalid full sampled parent matrix")
}
lookup <- readr::read_csv(private_result_path("network_case_map.csv"), show_col_types=FALSE)
require_columns(lookup, c("case_id", "id_int"), "private mapping")
if (nrow(lookup) != nrow(parents) || anyDuplicated(lookup$case_id) ||
    !identical(as.integer(lookup$id_int), seq_len(nrow(parents)))) stop("Incorrect network mapping")
match_idx <- match(as.character(d$case_id), as.character(lookup$case_id))
if (anyNA(match_idx)) stop("Unmatched regression case IDs; do NOT use row indices as a fallback")
# Check the outcome aligns to the exact 10k network mapping, not just correlation.
all_mean_ri <- tabulate(as.integer(parents[parents > 0L]), nbins=nrow(parents)) / ncol(parents)
max_delta <- max(abs(d$mean_Ri - all_mean_ri[match_idx]))
if (max_delta > 1e-5) stop("Stored mean_Ri does not match full-network offspring counts; max error=", max_delta)
set.seed(20260708L)
selected <- sort(sample(seq_len(ncol(parents)), size=300L, replace=FALSE))
form <- primary_formulas()$M3
form_sample <- stats::update(form, Ri_sampled ~ .)
fits <- lapply(seq_along(selected), function(i) {
  col_id <- selected[i]
  tryCatch({
    sampled <- parents[, col_id]
    d$Ri_sampled <- tabulate(sampled[sampled > 0L], nbins=nrow(parents))[match_idx]
    model <- AER::tobit(form_sample, data=d, left=0)
    list(ok=TRUE, message="", coefficients=tobit_coef(model, paste0("network_", col_id)),
         joint=all_joint_tests(model, paste0("network_", col_id)))
  }, error=function(e) list(ok=FALSE, message=conditionMessage(e)))
})
convergence <- data.frame(network_column=selected, success=vapply(fits, `[[`, logical(1), "ok"),
                          message=vapply(fits, `[[`, character(1), "message"))
safe_write_csv(convergence, "sampled_network_convergence.csv")
good <- Filter(function(x) x$ok, fits)
if (length(good) != length(selected)) stop("Some sampled network models failed; see private convergence log")
coef_all <- do.call(rbind, lapply(good, `[[`, "coefficients"))
joint_all <- do.call(rbind, lapply(good, `[[`, "joint"))
safe_write_csv(coef_all, "sampled_network_all_coefficients.csv")
safe_write_csv(joint_all, "sampled_network_all_joint_tests.csv")
linear_terms <- c("case_typeLocal","z_BI","z_outflow","z_dens_park","z_log_hosp","z_log_trans")
coef_focus <- subset(coef_all, term %in% linear_terms)
# Preserve the final uncertainty-table fields and percentage units (0..100).
linear_labels <- c(
  "case_typeLocal"="Locally acquired case", "z_BI"="Breteau Index",
  "z_outflow"="Population outflow",
  "z_dens_park"="Park and green-space count within 500 m",
  "z_log_hosp"="Healthcare facility count within 500 m",
  "z_log_trans"="Distance to the nearest public transport facility")
linear <- do.call(rbind,lapply(split(coef_focus,coef_focus$term),function(z) {
  data.frame(variable=unname(linear_labels[z$term[1]]),term=z$term[1],
    successful_models=nrow(z), median_coefficient=median(z$estimate),
    percentile_2_5=unname(quantile(z$estimate,.025)),
    percentile_97_5=unname(quantile(z$estimate,.975)),
    positive_direction_percent=100*mean(z$estimate>0),
    negative_direction_percent=100*mean(z$estimate<0),
    p_less_0_05_percent=100*mean(z$p_value<.05))
}))
spline <- do.call(rbind,lapply(split(joint_all,joint_all$term),function(z) {
  data.frame(variable=z$term[1],successful_models=nrow(z),
    median_wald_chisq=median(z$chi_square),
    percentile_2_5_wald=unname(quantile(z$chi_square,.025)),
    percentile_97_5_wald=unname(quantile(z$chi_square,.975)),
    median_p_value=median(z$p_value),
    p_less_0_05_percent=100*mean(z$p_value<.05))
}))
safe_write_csv(linear,"sampled_network_linear_stability.csv")
safe_write_csv(spline,"sampled_network_spline_stability.csv")
message("300 sampled-network refits complete; original percentage/quantile fields saved privately.")
