# Scenario 4: clean oracle verification of the boundary part of Theorem 3
#
# This script deliberately has one job.  It verifies the singular boundary
# rate in Theorem 3(ii) with K0 = 2 and one component at alpha = 1/2.  The
# mixing weight, the other location, scales, shapes and skew parameters are
# fixed at their population values.  Only the boundary location is optimized.
# Consequently the output is a theorem diagnostic, not a claim about the
# finite-sample behaviour of the unknown-nuisance/profile estimator.

.onf_script_file <- local({
  command <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", command, value = TRUE)
  candidate <- if (length(hit)) sub("^--file=", "", hit[1L]) else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(candidate) || !nzchar(candidate))
    candidate <- file.path(getwd(), "scenario4_boundary_oracle.R")
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
})
.onf_script_dir <- dirname(.onf_script_file)
.onf_engine <- file.path(.onf_script_dir, "scenario4_rates.R")
if (!file.exists(.onf_engine))
  stop("Required adjacent engine not found: ", .onf_engine)
source(.onf_engine, chdir = TRUE)

.onf_version <- "scenario4_boundary_oracle_v1"
.onf_design <- list(
  label = "Boundary oracle: K0=2, alpha1=1/2",
  seed = 43001L,
  pi = c(0.5, 0.5),
  mu = c(-3, 3),
  sigma = c(1, 1),
  alpha = c(0.5, 2),
  tau = c(0.5, 0.5),
  beta = 0.5,
  boundary_component = 1L,
  target = "Theorem 3(ii) boundary component"
)
.onf_truth <- .onf_design[c("pi", "mu", "sigma", "alpha", "tau")]
.onf_boundary_component <- as.integer(.onf_design$boundary_component)
.onf_boundary_C <- function(pars, k) {
  # Cbar for the alpha_k = 1/2 boundary component.  This is the constant in
  # Var{sqrt(n log n)(mu_hat_k-mu_k*)} = 1/Cbar.
  f0 <- mixture_aepd_grid(
    pars$mu[k], pars$pi, pars$mu, pars$sigma, pars$alpha, pars$tau
  )
  denom <- 32 * pars$sigma[k]^3 * f0
  value <- pars$pi[k]^2 *
    (pars$tau[k] * (1 - pars$tau[k]))^2 / denom
  as.numeric(value)
}
.onf_Cbar <- .onf_boundary_C(.onf_truth, .onf_boundary_component)
.onf_Vmu <- 1 / .onf_Cbar

.onf_message <- function(...)
  message(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)

.onf_empty <- function(b, n, reason, kde = NULL) list(
  b = b, n = n, ok = FALSE, reason = reason,
  mu_hat = NA_real_, error = NA_real_, z = NA_real_,
  objective_hat = NA_real_, objective_truth = NA_real_,
  interval_boundary = NA,
  grid_n = if (is.null(kde)) NA_integer_ else kde$grid_n,
  points_per_bandwidth = if (is.null(kde)) NA_real_ else
    kde$points_per_bandwidth,
  grid_capped = if (is.null(kde)) NA else kde$grid_capped,
  outside_fraction = if (is.null(kde)) NA_real_ else kde$outside_fraction,
  truth_grid_mass = NA_real_
)

.onf_one_replication <- function(b, n, config) {
  truth <- .onf_truth
  k <- .onf_boundary_component
  set.seed(config$seed_base + b)
  x <- rmix_aepd(
    n, truth$pi, truth$mu, truth$sigma, truth$alpha, truth$tau
  )
  kde <- tryCatch(make_kde_theory(
    x, .onf_design$beta, config$c_h, config$tail_p,
    config$per_bandwidth, config$grid_min, config$grid_max
  ), error = function(e) NULL)
  if (is.null(kde)) return(.onf_empty(b, n, "theory KDE failed"))
  if (isTRUE(kde$grid_capped))
    return(.onf_empty(b, n, "grid cap reached", kde))

  other <- setdiff(seq_along(truth$pi), k)
  fixed_density <- vapply(other, function(j) daepd(
    kde$grid, truth$mu[j], truth$sigma[j], truth$alpha[j], truth$tau[j]
  ), numeric(length(kde$grid)))
  if (length(other) == 1L) fixed_density <- matrix(fixed_density, ncol = 1L)
  fixed_part <- as.vector(fixed_density %*% truth$pi[other])
  objective <- function(mu) {
    boundary_density <- daepd(
      kde$grid, mu, truth$sigma[k], truth$alpha[k], truth$tau[k]
    )
    mixture <- truth$pi[k] * boundary_density + fixed_part
    2 - 2 * sum(kde$w * kde$sqrt_phat * sqrt(pmax(mixture, 0)))
  }

  gap <- min(abs(truth$mu[k] - truth$mu[other]))
  theoretical_se <- sqrt(.onf_Vmu / (n * log(n)))
  profile_radius <- min(0.40 * gap, config$profile_radius_se * theoretical_se)
  interval <- truth$mu[k] + c(-profile_radius, profile_radius)
  optimization_tolerance <- max(1e-10, config$tol_fraction * theoretical_se)
  scan_mu <- seq(interval[1L], interval[2L],
                 length.out = config$profile_grid_points)
  scan_objective <- vapply(scan_mu, objective, numeric(1))
  local_index <- which(
    scan_objective[2L:(length(scan_objective) - 1L)] <=
      scan_objective[1L:(length(scan_objective) - 2L)] &
    scan_objective[2L:(length(scan_objective) - 1L)] <=
      scan_objective[3L:length(scan_objective)]
  ) + 1L
  if (!length(local_index)) local_index <- which.min(scan_objective)
  refinements <- lapply(local_index, function(j) tryCatch(
    stats::optimize(
      objective, interval = scan_mu[c(j - 1L, j + 1L)],
      tol = optimization_tolerance
    ), error = function(e) NULL
  ))
  refinements <- Filter(Negate(is.null), refinements)
  candidates <- c(
    lapply(refinements, function(o)
      c(mu = o$minimum, objective = o$objective)),
    list(c(mu = truth$mu[k], objective = objective(truth$mu[k])))
  )
  values <- vapply(candidates, function(o) o[["objective"]], numeric(1))
  optimum <- candidates[[which.min(values)]]
  if (!all(is.finite(optimum)))
    return(.onf_empty(b, n, "one-dimensional optimization failed", kde))

  mu_hat <- optimum[["mu"]]
  objective_truth <- objective(truth$mu[k])
  at_boundary <- min(abs(mu_hat - interval)) <=
    10 * optimization_tolerance
  if (isTRUE(at_boundary))
    return(.onf_empty(b, n, "oracle minimum reached interval boundary", kde))

  truth_G <- component_matrix(kde$grid, truth)
  error <- mu_hat - truth$mu[k]
  list(
    b = b, n = n, ok = TRUE, reason = "success",
    mu_hat = mu_hat, error = error,
    z = error / theoretical_se,
    theoretical_se = theoretical_se,
    optimization_tolerance = optimization_tolerance,
    objective_hat = optimum[["objective"]],
    objective_truth = objective_truth,
    interval_boundary = FALSE,
    grid_n = kde$grid_n,
    points_per_bandwidth = kde$points_per_bandwidth,
    grid_capped = kde$grid_capped,
    outside_fraction = kde$outside_fraction,
    truth_grid_mass = sum(
      kde$w * as.vector(truth_G %*% truth$pi)
    )
  )
}

.onf_worker <- function(b)
  .onf_one_replication(b, .onf_n, .onf_config)

.onf_finite <- function(x) x[is.finite(x)]
.onf_stat <- function(x, fun) {
  x <- .onf_finite(x)
  if (length(x)) fun(x) else NA_real_
}

.onf_qq_values <- function(z, level = 0.95) {
  z <- sort(.onf_finite(z))
  m <- length(z)
  if (m < 3L) return(NULL)
  i <- seq_len(m)
  alpha <- 1 - level
  probability <- stats::ppoints(m)
  list(
    theoretical = stats::qnorm(probability), observed = z,
    lower = stats::qnorm(stats::qbeta(
      alpha / 2, shape1 = i, shape2 = m + 1L - i
    )),
    upper = stats::qnorm(stats::qbeta(
      1 - alpha / 2, shape1 = i, shape2 = m + 1L - i
    ))
  )
}

.onf_boundary_diagnostic_scaling <- function(n, c_h) {
  sigma_k <- .onf_truth$sigma[.onf_boundary_component]
  s_one_step <- .onf_Cbar^(-1 / 2)
  log_n <- log(n)
  kappa1 <- 1 + (
    log(log_n) - 2 * log(abs(s_one_step) / sigma_k)
  ) / log_n
  kappa2 <- 1 - 2 * log(c_h / sigma_k) / log_n
  if (!is.finite(kappa1) || !is.finite(kappa2) ||
      kappa1 <= 0 || kappa2 <= 0) {
    stop(
      "Finite-n boundary diagnostic scaling is undefined: ",
      "n=", n, ", kappa1=", signif(kappa1, 6),
      ", kappa2=", signif(kappa2, 6)
    )
  }
  list(
    kappa1 = kappa1, kappa2 = kappa2,
    Cbar_ratio = kappa1^2 / kappa2,
    z_multiplier = kappa1 / sqrt(kappa2)
  )
}

.onf_summary <- function(result) {
  rows <- lapply(result$by_n, function(record) {
    ok <- vapply(record$reps, function(r) isTRUE(r$ok), logical(1))
    error <- vapply(record$reps, function(r) r$error, numeric(1))
    z_raw <- vapply(record$reps, function(r) r$z, numeric(1))
    keep <- is.finite(error) & is.finite(z_raw)
    error <- error[keep]; z_raw <- z_raw[keep]
    n <- record$n
    c_h <- if (!is.null(record$config$c_h)) record$config$c_h else .8
    diagnostic <- .onf_boundary_diagnostic_scaling(n, c_h)
    z_diag <- z_raw * diagnostic$z_multiplier
    qq <- .onf_qq_values(z_diag)
    variance <- if (length(error) >= 2L) stats::var(error) else NA_real_
    mse <- if (length(error)) mean(error^2) else NA_real_
    skewness <- if (length(z_diag) >= 3L && stats::sd(z_diag) > 0) {
      centred <- z_diag - mean(z_diag)
      mean(centred^3) / stats::sd(z_diag)^3
    } else NA_real_
    excess_kurtosis <- if (length(z_diag) >= 4L && stats::var(z_diag) > 0) {
      centred <- z_diag - mean(z_diag)
      mean(centred^4) / stats::var(z_diag)^2 - 3
    } else NA_real_
    data.frame(
      n = n, B = record$B, valid = sum(ok & keep),
      bias = .onf_stat(error, mean), mse = mse, variance = variance,
      n_mse = n * mse, nlogn_mse = n * log(n) * mse,
      nlogn_variance = n * log(n) * variance,
      theory_variance = .onf_Vmu,
      Cbar_ratio = diagnostic$Cbar_ratio,
      diagnostic_variance = .onf_Vmu / diagnostic$Cbar_ratio,
      asymptotic_normalized_mse = n * log(n) * mse / .onf_Vmu,
      asymptotic_normalized_variance = n * log(n) * variance / .onf_Vmu,
      diagnostic_normalized_mse = n * log(n) * mse /
        (.onf_Vmu / diagnostic$Cbar_ratio),
      diagnostic_normalized_variance = n * log(n) * variance /
        (.onf_Vmu / diagnostic$Cbar_ratio),
      z_asym_mean = .onf_stat(z_raw, mean),
      z_asym_sd = if (length(z_raw) >= 2L) stats::sd(z_raw) else NA_real_,
      z_diag_mean = .onf_stat(z_diag, mean),
      z_diag_sd = if (length(z_diag) >= 2L) stats::sd(z_diag) else NA_real_,
      z_skewness = skewness, z_excess_kurtosis = excess_kurtosis,
      coverage_asymptotic = if (length(z_raw))
        mean(abs(z_raw) <= stats::qnorm(.975)) else NA_real_,
      coverage_diagnostic = if (length(z_diag))
        mean(abs(z_diag) <= stats::qnorm(.975)) else NA_real_,
      qq_outside_pointwise_band = if (is.null(qq)) NA_real_ else
        mean(qq$observed < qq$lower | qq$observed > qq$upper),
      grid_n_median = .onf_stat(vapply(
        record$reps, function(r) as.numeric(r$grid_n), numeric(1)
      ), stats::median),
      grid_n_max = .onf_stat(vapply(
        record$reps, function(r) as.numeric(r$grid_n), numeric(1)
      ), max),
      grid_capped = sum(vapply(
        record$reps, function(r) isTRUE(r$grid_capped), logical(1)
      )),
      outside_fraction_max = .onf_stat(vapply(
        record$reps, function(r) r$outside_fraction, numeric(1)
      ), max),
      truth_grid_mass_median = .onf_stat(vapply(
        record$reps, function(r) r$truth_grid_mass, numeric(1)
      ), stats::median),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.onf_plot_qq <- function(result, file, qq_n = NULL) {
  n_values <- sort(as.integer(names(result$by_n)))
  if (is.null(qq_n)) qq_n <- max(n_values)
  qq_n <- as.integer(qq_n)
  if (!qq_n %in% n_values)
    stop("qq_n is not present in result$by_n: ", qq_n)
  record <- result$by_n[[as.character(qq_n)]]
  z <- vapply(record$reps, function(r) r$z, numeric(1))
  diagnostic <- .onf_boundary_diagnostic_scaling(qq_n, record$config$c_h)
  z <- z * diagnostic$z_multiplier
  qq <- .onf_qq_values(z)
  grDevices::pdf(file, width = 6.4, height = 5.0)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (is.null(qq)) {
    graphics::plot.new()
    return(invisible(file))
  }
  ylim <- range(c(qq$observed, qq$lower, qq$upper), finite = TRUE)
  graphics::plot(
    qq$theoretical, qq$observed, type = "n", ylim = ylim,
    xlab = "Theoretical normal quantiles",
    ylab = "Finite-n corrected quantiles"
  )
  graphics::polygon(
    c(qq$theoretical, rev(qq$theoretical)),
    c(qq$lower, rev(qq$upper)), border = NA,
    col = grDevices::adjustcolor("grey70", alpha.f = 0.55)
  )
  graphics::abline(0, 1, col = "firebrick", lwd = 1.8)
  graphics::points(qq$theoretical, qq$observed, pch = 16, cex = 0.48)
  graphics::legend(
    "topleft", bty = "n",
    legend = c("95% pointwise band", "N(0,1) line", "estimated quantiles"),
    fill = c(grDevices::adjustcolor("grey70", alpha.f = 0.55), NA, NA),
    border = NA, lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
    pch = c(NA, NA, 16), col = c(NA, "firebrick", "black"), cex = .8
  )
  invisible(file)
}

.onf_plot_rates <- function(summary, file) {
  grDevices::pdf(file, width = 6.4, height = 5.0)
  on.exit(grDevices::dev.off(), add = TRUE)
  ylim <- range(c(1, summary$diagnostic_normalized_mse,
                  summary$diagnostic_normalized_variance),
                finite = TRUE)
  graphics::plot(
    summary$n, summary$diagnostic_normalized_mse, type = "b", pch = 16,
    lwd = 2,
    log = "x", ylim = ylim, col = "firebrick", xlab = "n",
    ylab = "Ratio to finite-n diagnostic variance"
  )
  graphics::lines(summary$n, summary$diagnostic_normalized_variance, type = "b",
                  pch = 17, lwd = 2, col = "steelblue")
  graphics::abline(h = 1, lty = 2, col = "grey35")
  graphics::legend(
    "topleft", bty = "n", lwd = 2, pch = c(16, 17),
    col = c("firebrick", "steelblue"),
    legend = c("n log(n) MSE / V[n]", "n log(n) Var / V[n]")
  )
  invisible(file)
}

.onf_write_outputs <- function(result, out_dir, qq_n = NULL) {
  summary <- .onf_summary(result)
  utils::write.csv(summary,
                   file.path(out_dir, "scenario4_boundary_oracle_summary.csv"),
                   row.names = FALSE)
  saveRDS(result,
          file.path(out_dir, "scenario4_boundary_oracle_results.rds"))
  .onf_plot_qq(result,
               file.path(out_dir, "scenario4_boundary_qq.pdf"), qq_n)
  .onf_plot_rates(
    summary, file.path(out_dir, "scenario4_boundary_rates.pdf")
  )
  summary
}

run_scenario4_boundary_oracle <- function(
    n_values = c(200L, 500L, 1000L, 2000L, 5000L, 10000L, 20000L),
    B = 500L, cores = 4L, c_h = 0.8,
    tail_p = 0, per_bandwidth = 4,
    grid_min = 512L, grid_max = 262144L,
    tol_fraction = 1e-4, profile_radius_se = 8,
    profile_grid_points = 65L,
    out_dir = file.path(.onf_script_dir, "scenario4_boundary_oracle_results"),
    cache_dir = file.path(.onf_script_dir, "oracle_nuisance_fixed_resultsB500"),
    qq_n = NULL, force = FALSE, quick = FALSE) {
  if (quick) {
    n_values <- c(200L, 500L); B <- 20L; cores <- 1L
    grid_max <- 32768L
  }
  if (!is.finite(tail_p) || tail_p < 0 || tail_p >= .5)
    stop("tail_p must lie in [0,0.5).")
  if (!is.finite(profile_radius_se) || profile_radius_se <= 0)
    stop("profile_radius_se must be positive.")
  if (!is.finite(profile_grid_points) || profile_grid_points < 5L)
    stop("profile_grid_points must be at least five.")
  profile_grid_points <- as.integer(profile_grid_points)
  if (profile_grid_points %% 2L == 0L)
    profile_grid_points <- profile_grid_points + 1L
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result <- list(
    version = .onf_version, design = .onf_design, truth = .onf_truth,
    Cbar = .onf_Cbar, Vmu = .onf_Vmu,
    n_values = as.integer(n_values), B = as.integer(B), by_n = list()
  )
  .onf_message(
    "Scenario 4 boundary check: K0=2; n=", paste(n_values, collapse = ","),
    "; B=", B, "; cores=", cores, "; tail_p=", tail_p,
    "; grid_max=", grid_max, "."
  )
  for (n_index in seq_along(n_values)) {
    n <- as.integer(n_values[n_index])
    cache_name <- sprintf(
      "boundary_oracle_n%d_B%d_ch%.2f_tail%.3g_grid%d_rad%.1f_scan%d.rds",
      n, B, c_h, tail_p, as.integer(grid_max), profile_radius_se,
      profile_grid_points)
    cache <- file.path(out_dir, cache_name)
    legacy_cache <- file.path(cache_dir, sprintf(
      "oracle_fixed_n%d_B%d_ch%.2f_tail%.3g_grid%d_rad%.1f_scan%d.rds",
      n, B, c_h, tail_p, as.integer(grid_max), profile_radius_se,
      profile_grid_points))
    if (file.exists(cache) && !force) {
      cached <- readRDS(cache)
      if (identical(cached$version, .onf_version)) {
        result$by_n[[as.character(n)]] <- cached
        .onf_message("[", n_index, "/", length(n_values), "] loaded ",
                     basename(cache))
        next
      }
    }
    if (!force && !file.exists(cache) && file.exists(legacy_cache)) {
      cached <- readRDS(legacy_cache)
      # The old cache was generated by the same one-dimensional oracle
      # estimator.  Reuse it only after checking the dimensions and config;
      # the new script changes labels and diagnostics, not the estimator.
      same_n <- identical(as.integer(cached$n), n)
      same_B <- identical(as.integer(cached$B), as.integer(B))
      same_config <- is.list(cached$config) &&
        isTRUE(all.equal(as.numeric(cached$config$c_h), as.numeric(c_h))) &&
        isTRUE(all.equal(as.numeric(cached$config$tail_p), as.numeric(tail_p)))
      if (same_n && same_B && same_config && is.list(cached$reps)) {
        cached$version <- .onf_version
        saveRDS(cached, cache)
        result$by_n[[as.character(n)]] <- cached
        .onf_message("[", n_index, "/", length(n_values), "] reused legacy ",
                     basename(legacy_cache))
        next
      }
    }
    config <- list(
      seed_base = .onf_design$seed + 200000L * n_index,
      c_h = c_h, tail_p = tail_p, per_bandwidth = per_bandwidth,
      grid_min = as.integer(grid_min), grid_max = as.integer(grid_max),
      tol_fraction = tol_fraction,
      profile_radius_se = profile_radius_se,
      profile_grid_points = as.integer(profile_grid_points)
    )
    assign(".onf_n", n, envir = globalenv())
    assign(".onf_config", config, envir = globalenv())
    started <- proc.time()[["elapsed"]]
    report_every <- max(1L, ceiling(B / 10L))
    report <- function(done, reps) {
      valid <- sum(vapply(reps, function(r) isTRUE(r$ok), logical(1)))
      elapsed <- proc.time()[["elapsed"]] - started
      eta <- elapsed * (B - done) / done
      .onf_message(sprintf(
        "[%d/%d] n=%d: %d/%d, valid=%d, elapsed=%.1f min, ETA=%.1f min.",
        n_index, length(n_values), n, done, B, valid,
        elapsed / 60, eta / 60
      ))
    }
    if (cores > 1L) {
      cl <- parallel::makeCluster(min(as.integer(cores), B))
      reps <- tryCatch({
        parallel::clusterExport(
          cl, c(".onf_script_file", ".onf_n", ".onf_config"),
          envir = globalenv()
        )
        parallel::clusterEvalQ(cl, source(.onf_script_file, chdir = TRUE))
        output <- vector("list", B)
        chunks <- split(seq_len(B), ceiling(seq_len(B) / report_every))
        for (chunk in chunks) {
          output[chunk] <- parallel::parLapplyLB(cl, chunk, .onf_worker)
          report(max(chunk), output[seq_len(max(chunk))])
        }
        output
      }, finally = try(parallel::stopCluster(cl), silent = TRUE))
    } else {
      reps <- vector("list", B)
      for (b in seq_len(B)) {
        reps[[b]] <- .onf_worker(b)
        if (b == 1L || b %% report_every == 0L || b == B)
          report(b, reps[seq_len(b)])
      }
    }
    record <- list(
      version = .onf_version, n = n, B = B, config = config,
      reps = reps, elapsed = proc.time()[["elapsed"]] - started
    )
    saveRDS(record, cache)
    result$by_n[[as.character(n)]] <- record
  }
  .onf_message("Writing summary and PDF diagnostics.")
  if (is.null(qq_n)) qq_n <- max(as.integer(n_values))
  result$qq_n <- as.integer(qq_n)
  summary <- .onf_write_outputs(result, out_dir, qq_n)
  .onf_message("Completed: ", out_dir)
  invisible(c(result, list(summary = summary)))
}

.onf_parse_cli <- function(arguments) {
  config <- list(
    n_values = c(200L, 500L, 1000L, 2000L, 5000L, 10000L, 20000L),
    B = 500L, cores = 4L, c_h = .8, tail_p = 0,
    grid_max = 262144L, tol_fraction = 1e-4,
    profile_radius_se = 8, profile_grid_points = 65L,
    out_dir = file.path(.onf_script_dir, "scenario4_boundary_oracle_results"),
    cache_dir = file.path(.onf_script_dir, "oracle_nuisance_fixed_resultsB500"),
    qq_n = NULL, force = FALSE, quick = FALSE
  )
  value <- function(a, prefix) sub(paste0("^", prefix), "", a)
  for (a in arguments) {
    if (a == "--quick") config$quick <- TRUE
    else if (a == "--force") config$force <- TRUE
    else if (grepl("^--n=", a))
      config$n_values <- as.integer(strsplit(value(a, "--n="), ",")[[1L]])
    else if (grepl("^--B=", a)) config$B <- as.integer(value(a, "--B="))
    else if (grepl("^--cores=", a))
      config$cores <- as.integer(value(a, "--cores="))
    else if (grepl("^--c-h=", a))
      config$c_h <- as.numeric(value(a, "--c-h="))
    else if (grepl("^--tail-p=", a))
      config$tail_p <- as.numeric(value(a, "--tail-p="))
    else if (grepl("^--grid-max=", a))
      config$grid_max <- as.integer(value(a, "--grid-max="))
    else if (grepl("^--tol-fraction=", a))
      config$tol_fraction <- as.numeric(value(a, "--tol-fraction="))
    else if (grepl("^--profile-radius-se=", a))
      config$profile_radius_se <- as.numeric(value(a, "--profile-radius-se="))
    else if (grepl("^--profile-grid-points=", a))
      config$profile_grid_points <- as.integer(
        value(a, "--profile-grid-points=")
      )
    else if (grepl("^--out-dir=", a))
      config$out_dir <- value(a, "--out-dir=")
    else if (grepl("^--cache-dir=", a))
      config$cache_dir <- value(a, "--cache-dir=")
    else if (grepl("^--qq-n=", a))
      config$qq_n <- as.integer(value(a, "--qq-n="))
    else stop("Unknown option: ", a)
  }
  config
}

if (sys.nframe() == 0L) {
  .onf_result <- do.call(
    run_scenario4_boundary_oracle,
    .onf_parse_cli(commandArgs(trailingOnly = TRUE))
  )
  print(.onf_result$summary, row.names = FALSE)
}
