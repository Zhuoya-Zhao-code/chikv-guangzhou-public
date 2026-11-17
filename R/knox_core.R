# Knox permutation test. Compute all thresholds on the SAME permutations:
# protects against inconsistent Monte Carlo draws across cells.
knox_grid <- function(xy_m, onset_dates,
                      spatial_m = c(50, 100, 200, 300, 500, 1000, 1500),
                      temporal_days = c(7, 14, 21, 28),
                      permutations = 999L, seed = 123L) {
  xy_m <- as.matrix(xy_m)
  stopifnot(nrow(xy_m) == length(onset_dates), ncol(xy_m) == 2L,
            permutations >= 1, length(onset_dates) >= 3)
  onset <- as.numeric(onset_dates)
  if (any(!is.finite(onset)) || any(!is.finite(xy_m))) stop("Invalid coordinates or dates")
  ds <- as.numeric(stats::dist(xy_m))
  dt <- as.numeric(stats::dist(onset))
  grid <- expand.grid(spatial_m = spatial_m, temporal_days = temporal_days)
  spatial_idx <- lapply(spatial_m, function(z) which(ds <= z))
  observed <- mapply(function(s, t) sum(dt[spatial_idx[[match(s, spatial_m)]]] <= t),
                     grid$spatial_m, grid$temporal_days)
  perm_counts <- matrix(0L, nrow = nrow(grid), ncol = permutations)
  set.seed(seed)
  for (p in seq_len(permutations)) {
    p_dt <- as.numeric(stats::dist(sample(onset)))
    perm_counts[, p] <- mapply(function(s, t) sum(p_dt[spatial_idx[[match(s, spatial_m)]]] <= t),
                                grid$spatial_m, grid$temporal_days)
  }
  expected <- rowMeans(perm_counts)
  grid$observed_pairs <- observed
  grid$expected_pairs <- expected
  grid$observed_expected_ratio <- ifelse(expected > 0, observed / expected, NA_real_)
  grid$permutation_p <- (rowSums(perm_counts >= observed) + 1) / (permutations + 1)
  grid
}
