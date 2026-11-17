# Final onset-/residential-distance-compatible reconstruction; no intervention penalty.
# For the primary 10,000-network run, preserve the original simulation-major
# RNG loop order from "1 传播链重建.Rmd". Sensitivity scenarios use the distinct,
# case-major order documented in the late sensitivity block of that Rmd.
network_parameters <- function() list(
  si_mean = 8, si_sd = 3, distance_scale_m = 250,
  distance_cutoff_m = 1300, max_lag_days = 30L,
  background_weight = 1e-4, n_samples = 10000L, seed = 2026L
)
validate_network_input <- function(dat) {
  require_columns(dat, c("case_id", "onset_date", "case_type", "utm_x", "utm_y"), "network cases")
  check_case_ids(dat$case_id)
  dat$onset_date <- make_onset_date(dat$onset_date)
  dat$case_type <- tolower(trimws(as.character(dat$case_type)))
  if (anyNA(dat$case_type) || !all(dat$case_type %in% c("local", "imported"))) {
    stop("case_type must be Local or Imported")
  }
  dat$utm_x <- suppressWarnings(as.numeric(dat$utm_x))
  dat$utm_y <- suppressWarnings(as.numeric(dat$utm_y))
  if (any(!is.finite(dat$utm_x)) || any(!is.finite(dat$utm_y))) {
    stop("Coordinates must be finite UTM metres")
  }
  # The source sorted by onset date; ties must preserve source-input order.
  dat <- dat[order(dat$onset_date, seq_len(nrow(dat))), , drop = FALSE]
  dat$id_int <- seq_len(nrow(dat))
  dat
}
make_score_matrix <- function(dat, params = network_parameters()) {
  if (params$si_mean <= 0 || params$si_sd <= 0 || params$distance_scale_m <= 0 ||
      params$distance_cutoff_m <= 0 || params$max_lag_days < 1 ||
      !is.finite(params$background_weight) || params$background_weight <= 0) {
    stop("Invalid reconstruction parameters")
  }
  n <- nrow(dat)
  dates <- dat$onset_date
  xy <- as.matrix(dat[, c("utm_x", "utm_y")])
  shape <- (params$si_mean / params$si_sd)^2
  scale <- params$si_sd^2 / params$si_mean
  scores <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (dat$case_type[i] != "local") next
    dt <- as.numeric(dates[i] - dates)
    candidates <- which(dt > 0 & dt <= params$max_lag_days)
    if (!length(candidates)) next
    distance <- sqrt((xy[i, 1] - xy[candidates, 1])^2 +
                     (xy[i, 2] - xy[candidates, 2])^2)
    eligible <- distance <= params$distance_cutoff_m
    if (!any(eligible)) next
    upstream <- candidates[eligible]
    scores[i, upstream] <- stats::dgamma(dt[upstream], shape = shape, scale = scale) *
      exp(-distance[eligible] / params$distance_scale_m)
  }
  scores
}
sample_networks <- function(dat, scores, params = network_parameters(),
                            rng_order = c("simulation_major", "case_major")) {
  rng_order <- match.arg(rng_order)
  n <- nrow(dat)
  if (!identical(dim(scores), c(n, n))) stop("Score-matrix dimensions do not match cases")
  if (anyNA(scores) || any(!is.finite(scores)) || any(scores < 0)) stop("Invalid score matrix")
  nsim <- as.integer(params$n_samples)
  if (is.na(nsim) || nsim < 1L) stop("n_samples must be a positive integer")
  parents <- matrix(0L, nrow = n, ncol = nsim)
  ids <- c(seq_len(n), 0L)
  local <- which(dat$case_type == "local")
  weights <- lapply(local, function(i) c(scores[i, ], params$background_weight))
  set.seed(params$seed)
  if (rng_order == "simulation_major") {
    # Primary analysis: for(s in 1:n_sims) for(i in 1:n_cases).
    for (s in seq_len(nsim)) {
      for (k in seq_along(local)) {
        parents[local[k], s] <- sample(ids, size = 1L, prob = weights[[k]])
      }
    }
  } else {
    # Late sensitivity section: for(i in local_indices) sample(..., size=n_sims,
    # replace=TRUE). Do not silently change this order: it changes the draws.
    for (k in seq_along(local)) {
      parents[local[k], ] <- sample(ids, size = nsim, replace = TRUE, prob = weights[[k]])
    }
  }
  parents
}
# A per-case modal assignment over ALL states (including unknown/background=0).
# This is the late Rmd's descriptive maximum-marginal-posterior map, NOT a
# maximization over the joint probability of one entire network.
# Ties choose 0 where present, otherwise the smallest numeric source index.
posterior_source_summary <- function(parents) {
  if (!is.matrix(parents) || !nrow(parents) || !ncol(parents) || anyNA(parents)) {
    stop("Full nonempty parent matrix required")
  }
  n <- nrow(parents)
  if (any(parents < 0 | parents > n | parents != as.integer(parents))) {
    stop("Illegal parent indices")
  }
  row_list <- lapply(seq_len(n), function(i) {
    x <- as.integer(parents[i, ])
    counts <- tabulate(x + 1L, nbins = n + 1L)
    largest <- max(counts)
    ties <- which(counts == largest) - 1L
    c(parent = min(ties), support_prob = largest / length(x),
      n_tied_sources = length(ties), no_observed_upstream_prob = mean(x == 0L))
  })
  result <- as.data.frame(do.call(rbind, row_list))
  result$parent <- as.integer(result$parent)
  result$n_tied_sources <- as.integer(result$n_tied_sources)
  result
}
modal_parent_with_background <- function(parents) posterior_source_summary(parents)$parent
individual_ri <- function(parents) {
  n <- nrow(parents)
  nsim <- ncol(parents)
  ri <- matrix(0L, nrow = n, ncol = nsim)
  for (s in seq_len(nsim)) {
    source <- parents[, s]
    ri[, s] <- tabulate(source[source > 0L], nbins = n)
  }
  ri
}
extract_network_summary <- function(dat, parents, max_distance_m = 1300) {
  if (nrow(dat) != nrow(parents)) stop("Case-to-network shape mismatch")
  n <- nrow(dat)
  posterior <- posterior_source_summary(parents)
  selected <- posterior$parent
  ri <- individual_ri(parents)
  x <- dat[, c("case_id", "id_int", "case_type", "onset_date", "utm_x", "utm_y"), drop = FALSE]
  x$top_parent_id <- selected
  x$support_prob <- posterior$support_prob
  x$n_tied_sources <- posterior$n_tied_sources
  x$no_observed_upstream_prob <- posterior$no_observed_upstream_prob
  x$isolation_prob <- posterior$no_observed_upstream_prob
  x$mean_Ri <- rowMeans(ri)
  x$Ri_LCI <- apply(ri, 1, stats::quantile, probs = 0.025, names = FALSE)
  x$Ri_UCI <- apply(ri, 1, stats::quantile, probs = 0.975, names = FALSE)
  x$infect_distance_m <- NA_real_
  x$infect_time_diff_d <- NA_real_
  kids <- which(selected > 0L)
  if (length(kids)) {
    adults <- selected[kids]
    x$infect_distance_m[kids] <- sqrt((x$utm_x[kids] - x$utm_x[adults])^2 +
                                     (x$utm_y[kids] - x$utm_y[adults])^2)
    x$infect_time_diff_d[kids] <- as.numeric(x$onset_date[kids] - x$onset_date[adults])
    if (any(x$infect_time_diff_d[kids] <= 0) ||
        any(x$infect_distance_m[kids] > max_distance_m + 1e-8)) {
      stop("Chosen edge violates eligibility criteria")
    }
  }
  x$downstream_count <- tabulate(selected[selected > 0L], nbins = n)
  x$has_observed_upstream <- selected != 0L
  x$has_observed_downstream <- x$downstream_count > 0L
  x$is_unconnected <- !x$has_observed_upstream & !x$has_observed_downstream
  x$network_role <- ifelse(!x$has_observed_upstream & x$has_observed_downstream, "Seed",
                           ifelse(x$has_observed_upstream & x$has_observed_downstream, "Intermediate",
                                  ifelse(x$has_observed_upstream, "Terminal", "Unconnected")))
  if (sum(x$downstream_count) != length(kids) ||
      abs(sum(x$mean_Ri) - sum(parents > 0L) / ncol(parents)) > 1e-7) {
    stop("Network offspring accounting failed")
  }
  x
}
