# Weak components of a directed parent forest are exactly its ancestral-root
# groups. All singleton nodes remain in membership; only size >= 2 are clusters.
parent_root_depth <- function(parents) {
  n <- length(parents)
  if (anyNA(parents) || any(parents < 0 | parents > n | parents != as.integer(parents))) {
    stop("Invalid parent indices")
  }
  roots <- integer(n)
  depth <- integer(n)
  for (i in seq_len(n)) {
    at <- i
    steps <- 0L
    seen <- logical(n)
    repeat {
      if (seen[at]) stop("Cycle in descriptive parent network")
      seen[at] <- TRUE
      p <- as.integer(parents[at])
      if (p == 0L) { roots[i] <- at; depth[i] <- steps; break }
      steps <- steps + 1L
      at <- p
    }
  }
  data.frame(seed_id = roots, generation_depth = depth)
}
find_cluster_root <- function(parents) parent_root_depth(parents)$seed_id
cluster_statistics <- function(dat) {
  require_columns(dat, c("case_id", "id_int", "top_parent_id", "onset_date", "utm_x", "utm_y"),
                  "network summary")
  dat <- dat[order(dat$id_int), , drop = FALSE]
  n <- nrow(dat)
  if (!identical(as.integer(dat$id_int), seq_len(n))) stop("id_int must be ordered 1..N")
  dat$onset_date <- make_onset_date(dat$onset_date)
  tracing <- parent_root_depth(as.integer(dat$top_parent_id))
  dat$seed_id <- tracing$seed_id
  dat$generation_depth <- tracing$generation_depth
  dat$downstream_count <- tabulate(dat$top_parent_id[dat$top_parent_id > 0L], nbins = n)
  dat$has_observed_upstream <- dat$top_parent_id != 0L
  dat$has_observed_downstream <- dat$downstream_count > 0L
  dat$is_unconnected <- !dat$has_observed_upstream & !dat$has_observed_downstream
  grouped <- split(seq_len(n), dat$seed_id)
  rows <- lapply(grouped, function(indices) {
    xy <- as.matrix(dat[indices, c("utm_x", "utm_y")])
    data.frame(seed_id = dat$seed_id[indices[1]], cluster_size = length(indices),
               start_onset_date = min(dat$onset_date[indices]),
               end_onset_date = max(dat$onset_date[indices]),
               duration_days = as.numeric(max(dat$onset_date[indices]) - min(dat$onset_date[indices])),
               residential_diameter_m = if (length(indices) > 1L) max(stats::dist(xy)) else 0,
               maximum_generation_depth = max(dat$generation_depth[indices]))
  })
  all_components <- do.call(rbind, rows)
  clusters <- subset(all_components, cluster_size >= 2L)
  clusters <- clusters[order(-clusters$cluster_size, clusters$seed_id), , drop=FALSE]
  rownames(clusters) <- NULL
  clusters$cluster_id <- seq_len(nrow(clusters))
  lookup <- setNames(clusters$cluster_id, as.character(clusters$seed_id))
  dat$cluster_id <- as.integer(lookup[as.character(dat$seed_id)])
  if (sum(dat$is_unconnected) != sum(all_components$cluster_size == 1L) ||
      sum(clusters$cluster_size) + sum(dat$is_unconnected) != n ||
      sum(dat$top_parent_id != 0L) != sum(clusters$cluster_size) - nrow(clusters)) {
    stop("Singleton/connected-component accounting mismatch")
  }
  list(membership = dat[, c("case_id", "id_int", "seed_id", "cluster_id", "generation_depth",
                             "downstream_count", "is_unconnected")],
       clusters = clusters, all_components = all_components)
}
# Mirrors the late Rmd OLS: all clusters of size >= 2, including zero-duration ones.
cluster_slope <- function(clusters) {
  if (nrow(clusters) < 3L || length(unique(clusters$duration_days)) < 2L) {
    return(data.frame(estimate=NA_real_, lower=NA_real_, upper=NA_real_,
                      p_value=NA_real_, r_squared=NA_real_))
  }
  fit <- stats::lm(residential_diameter_m ~ duration_days, data=clusters)
  ci <- stats::confint(fit, "duration_days", level=0.95)
  data.frame(estimate=unname(stats::coef(fit)["duration_days"]),
             lower=unname(ci[1]), upper=unname(ci[2]),
             p_value=unname(summary(fit)$coefficients["duration_days", "Pr(>|t|)"]),
             r_squared=summary(fit)$r.squared)
}
# Wilson score interval, replacing binom::binom.confint without a dependency.
wilson_ci <- function(x, n, z = stats::qnorm(.975)) {
  if (n <= 0 || x < 0 || x > n) return(c(lower=NA_real_, upper=NA_real_))
  p <- x/n
  denom <- 1+z^2/n
  centre <- (p+z^2/(2*n))/denom
  half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/denom
  c(lower=max(0,centre-half),upper=min(1,centre+half))
}
