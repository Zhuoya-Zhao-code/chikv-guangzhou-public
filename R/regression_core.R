# Current, temporally adjusted nested Tobit specification from the final scripts.
# Inputs are already derived / standardized in restricted data preparation.
prepare_regression <- function(dat) {
  required <- c("mean_Ri", "z_age", "sex", "occ_recode", "case_type", "Reporting_Delay",
                "matched_P", "z_BI", "z_outflow", "z_dens_park", "z_log_hosp", "z_log_trans",
                "log_dist_market", "utm_x_km", "utm_y_km", "onset_date")
  require_columns(dat, required, "Case_Data_Phase1_Ready.csv")
  if ("case_id" %in% names(dat)) check_case_ids(dat$case_id)
  dat$onset_date <- make_onset_date(dat$onset_date)
  dat$onset_day <- as.numeric(dat$onset_date - min(dat$onset_date))
  dat$sex <- factor(as.character(dat$sex), levels = c("Female", "Male"))
  dat$case_type <- factor(as.character(dat$case_type), levels = c("Imported", "Local"))
  dat$occ_recode <- factor(as.character(dat$occ_recode), levels = c(
    "Indoor/Sedentary", "Student/Childcare", "High Exposure/Mobile", "Home-based", "Other/Unknown"))
  covars <- c(required, "onset_day")
  if (any(!stats::complete.cases(dat[, covars]))) stop("Missing/unrecognized model values: cannot silently drop cases")
  if (any(!is.finite(dat$mean_Ri)) || any(dat$mean_Ri < 0)) stop("mean_Ri must be finite and >= 0")
  if (nrow(dat) != 1280L) warning("Expected 1280 complete regression cases; got ", nrow(dat),
                               ". Verify sample selection against the manuscript.")
  dat
}
primary_formulas <- function() {
  base <- c("z_age", "sex", "occ_recode", "case_type",
            "splines::ns(Reporting_Delay, df = 3)", "splines::ns(onset_day, df = 4)",
            "stats::poly(utm_x_km, utm_y_km, degree = 2)")
  env <- c("splines::ns(matched_P, df = 3)", "z_BI", "z_outflow")
  built <- c("z_dens_park", "z_log_hosp", "z_log_trans", "splines::ns(log_dist_market, df = 3)")
  fm <- function(terms) stats::as.formula(paste("mean_Ri ~", paste(terms, collapse = " + ")))
  list(M1 = fm(base), M2 = fm(c(base, env)), M3 = fm(c(base, env, built)))
}
tobit_coef <- function(model, model_name) {
  beta <- stats::coef(model)
  vc <- stats::vcov(model)
  se <- sqrt(diag(vc)[names(beta)])
  z <- beta / se
  data.frame(model = model_name, term = names(beta), estimate = unname(beta),
             std_error = unname(se), p_value = unname(2*stats::pnorm(abs(z), lower.tail = FALSE)),
             ci_lower = unname(beta - stats::qnorm(.975)*se),
             ci_upper = unname(beta + stats::qnorm(.975)*se))
}
wald_term <- function(model, label, match) {
  beta <- stats::coef(model); vc <- stats::vcov(model)
  hit <- grepl(match, names(beta), fixed = TRUE)
  if (!any(hit)) return(data.frame(term = label, df = NA_integer_, chi_square = NA_real_, p_value = NA_real_))
  x <- beta[hit]; v <- vc[names(x), names(x), drop = FALSE]
  df <- qr(v)$rank
  if (df != length(x)) stop("Rank-deficient coefficient covariance for ", label)
  chi <- as.numeric(t(x) %*% solve(v) %*% x)
  data.frame(term = label, df = df, chi_square = chi,
             p_value = stats::pchisq(chi, df, lower.tail = FALSE))
}
all_joint_tests <- function(model, model_name) {
  defs <- c("Reporting delay" = "Reporting_Delay", "Symptom onset date" = "onset_day",
            "Environmental Risk Index" = "matched_P", "Distance to the nearest market" = "log_dist_market")
  result <- do.call(rbind, Map(function(label, match) wald_term(model, label, match), names(defs), unname(defs)))
  result$model <- model_name
  result
}
model_fit_stats <- function(model, nul, model_name, n) {
  ll <- as.numeric(stats::logLik(model)); ll0 <- as.numeric(stats::logLik(nul))
  k <- attr(stats::logLik(model), "df")
  data.frame(model = model_name, n = n, log_likelihood = ll,
             AIC = -2*ll + 2*k, BIC = -2*ll + log(n)*k,
             McFadden_R2 = 1 - ll/ll0)
}
