# Primary nested models and selected sensitivity checks.
# Final onset-date-adjusted specifications, not the pre-adjustment Rmd models.
source("R/common.R")
source("R/regression_core.R")
if (!requireNamespace("AER", quietly = TRUE)) stop("Install R package AER")
if (!requireNamespace("spdep", quietly = TRUE)) stop("Install R package spdep")
d <- prepare_regression(read_private_csv("Case_Data_Phase1_Ready.csv"))
forms <- primary_formulas()
mods <- lapply(forms, function(f) AER::tobit(f, data = d, left = 0))
null <- AER::tobit(mean_Ri ~ 1, data = d, left = 0)
for (m in names(mods)) {
  if (length(stats::residuals(mods[[m]])) != nrow(d)) stop("Model silently omitted records: ", m)
}
metrics <- do.call(rbind, Map(function(m, name) model_fit_stats(m, null, name, nrow(d)), mods, names(mods)))
coeff <- do.call(rbind, Map(tobit_coef, mods, names(mods)))
joint <- do.call(rbind, Map(all_joint_tests, mods, names(mods)))
likelihood_ratio <- do.call(rbind, lapply(list(c("M1","M2"), c("M2","M3")), function(pair) {
  s <- mods[[pair[1]]]; l <- mods[[pair[2]]]
  df <- attr(stats::logLik(l), "df") - attr(stats::logLik(s), "df")
  stat <- 2*(as.numeric(stats::logLik(l))-as.numeric(stats::logLik(s)))
  if (df <= 0 || stat < -1e-8) stop("Invalid nested-model likelihood ratio")
  data.frame(comparison = paste(pair, collapse = " vs "), df = df,
             LR_chisq = max(stat, 0), p_value = stats::pchisq(max(stat, 0), df, lower.tail = FALSE))
}))
safe_write_csv(metrics, "primary_nested_model_metrics.csv")
safe_write_csv(coeff, "primary_coefficients.csv")
safe_write_csv(joint, "primary_joint_wald_tests.csv")
safe_write_csv(likelihood_ratio, "primary_likelihood_ratio_tests.csv")

# Right-boundary sensitivity: the last 30 days cannot have full onward follow-up.
cutoff <- as.Date("2025-11-30") - 30L
d_complete <- subset(d, onset_date <= cutoff)
if (nrow(d_complete) > 30) {
  complete_fit <- AER::tobit(forms$M3, data = d_complete, left = 0)
  safe_write_csv(tobit_coef(complete_fit, "M3_complete_30_day_window"), "sensitivity_complete_window_coefficients.csv")
  safe_write_csv(all_joint_tests(complete_fit, "M3_complete_30_day_window"), "sensitivity_complete_window_joint.csv")
}
# Local-only sensitivity: avoid a constant case_type term.
loc <- subset(d, case_type == "Local")
if (nrow(loc) > 30) {
  f_local <- stats::update(forms$M3, . ~ . - case_type)
  local_fit <- AER::tobit(f_local, data = loc, left = 0)
  safe_write_csv(tobit_coef(local_fit, "Local_only_M3"), "sensitivity_local_only_coefficients.csv")
  safe_write_csv(all_joint_tests(local_fit, "Local_only_M3"), "sensitivity_local_only_joint.csv")
}
# Diagnostic GVIF for the ACTUAL fitted Tobit model (not an LM proxy).
if (requireNamespace("car", quietly = TRUE)) {
  v <- tryCatch(car::vif(mods$M3), error = function(e) e)
  if (inherits(v, "error")) {
    warning("Tobit GVIF unavailable: ", conditionMessage(v))
  } else if (is.matrix(v)) {
    if (ncol(v) < 3L) stop("Unexpected car::vif matrix shape")
    safe_write_csv(data.frame(term=rownames(v), GVIF=v[,1], df=v[,2],
                              adjusted_GVIF = v[,3]), "diagnostic_gvif.csv")
  } else {
    safe_write_csv(data.frame(term=names(v), GVIF=as.numeric(v), df=1L,
                              adjusted_GVIF=sqrt(as.numeric(v))), "diagnostic_gvif.csv")
  }
}
# Model residual spatial autocorrelation; distances are kilometres.
resid <- as.numeric(stats::residuals(mods$M3, type = "response"))
coord <- as.matrix(d[, c("utm_x_km", "utm_y_km")])
moran_test <- function(kind, value, label, seed) {
  nb <- if (kind == "distance") {
    spdep::dnearneigh(coord, 0, value, longlat = FALSE)
  } else {
    spdep::knn2nb(spdep::knearneigh(coord, k = value, longlat = FALSE), sym = TRUE)
  }
  w <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  set.seed(seed)
  stat <- spdep::moran.mc(resid, w, nsim = 999L, zero.policy = TRUE)
  data.frame(weight_scheme = label, moran_I = as.numeric(stat$statistic),
             permutation_p = stat$p.value, isolated_cases = sum(spdep::card(nb)==0),
             disconnected_subgraphs = spdep::n.comp.nb(nb)$nc)
}
moran <- do.call(rbind, list(
  moran_test("distance", 1, "1000m", 124),
  moran_test("distance", 1.5, "1500m", 125),
  moran_test("distance", 2, "2000m", 126),
  moran_test("knn", 4, "symmetric 4NN", 127),
  moran_test("knn", 8, "symmetric 8NN", 128)
))
safe_write_csv(moran, "diagnostic_moran.csv")

# Figure 3B/C: marginally standardized LATENT Tobit predictor (not observed Ri).
# Change target variable for all individuals, average linear predictions;
# CI from the average design vector. No patient-level prediction is exported.
marginal_latent_curve <- function(model, data, target, values) {
  terms <- stats::delete.response(stats::terms(model))
  beta <- stats::coef(model)
  V <- stats::vcov(model)[names(beta), names(beta), drop = FALSE]
  do.call(rbind, lapply(values, function(v) {
    dd <- data; dd[[target]] <- v
    # Preserve fitting-time spline knots and polynomial coefficients by using
    # the fitted terms/predvars; never recalculate bases on a fixed target.
    mf <- stats::model.frame(terms, data = dd, na.action = stats::na.pass,
                             xlev = model$xlevels)
    X <- stats::model.matrix(terms, data = mf, contrasts.arg = model$contrasts)
    if (!all(names(beta) %in% colnames(X))) stop("Prediction matrix lacks fitted columns")
    X <- X[, names(beta), drop = FALSE]
    # Fail closed if the design construction differs from predict.survreg:
    # otherwise a visually plausible curve could have the wrong spline basis.
    direct <- as.numeric(stats::predict(model, newdata = dd, type = "lp"))
    reconstructed <- as.numeric(X %*% beta)
    if (length(direct) != nrow(X) || any(!is.finite(direct)) ||
        max(abs(direct - reconstructed)) > 1e-5) {
      stop("Prediction-design mismatch; check cached spline/poly coefficients")
    }
    a <- colMeans(X)
    mu <- sum(a * beta)
    se <- sqrt(as.numeric(t(a) %*% V %*% a))
    data.frame(variable=target, x=v, latent_mean=mu,
               lower=mu-1.959963985*se, upper=mu+1.959963985*se)
  }))
}
for (var in c("onset_day", "log_dist_market")) {
  values <- seq(stats::quantile(d[[var]], .01), stats::quantile(d[[var]], .99), length.out=100L)
  safe_write_csv(marginal_latent_curve(mods$M3, d, var, values), paste0("figure3_", var, "_latent_curve.csv"))
}
# Save only model objects in git-ignored private output.
saveRDS(mods$M3, private_result_path("model3_PRIVATE.rds"))
message("Primary and sensitivity Tobit analyses done; compare outputs with manuscript.")
