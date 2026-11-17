# Optional robustness: the source's later Tweedie sensitivity model.
# It is NOT a replacement for the primary Tobit analysis; coefficients use
# a LOG link and must NOT be compared directly to latent Tobit betas.
source("R/common.R")
source("R/regression_core.R")
if (!requireNamespace("mgcv",quietly=TRUE)) stop("Install mgcv")
d <- prepare_regression(read_private_csv("Case_Data_Phase1_Ready.csv"))
form <- primary_formulas()$M3
mod <- mgcv::gam(formula=form, data=d, family=mgcv::tw(link="log"),method="REML")
ct <- as.data.frame(summary(mod)$p.table)
names(ct) <- c("estimate","standard_error","statistic","p_value")
ct$term <- rownames(ct)
ct$mean_ratio <- exp(ct$estimate)
safe_write_csv(ct, "optional_tweedie_coefficients.csv")
results <- all_joint_tests(mod, "Tweedie_M3")
safe_write_csv(results, "optional_tweedie_joint.csv")
message("Optional Tweedie analysis done. Interpret as mean ratios under log link.")
