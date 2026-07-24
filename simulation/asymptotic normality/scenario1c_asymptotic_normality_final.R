# scenario4_asymptotic_normality_final.R
# ============================================================================
# Finite-sample check of the Gaussian approximation of Theorem 3(i) for the
# PMHD-MCP-R active-submodel minimizer. This is the experiment behind the
# asymptotic-normality / coverage table in the paper.
#
# SELF-CONTAINED: this script depends only on the core algorithm library
# simulation_stability_elbow_parallel_main.R (for daepd, rmix_aepd,
# pmhd_mcp_fit_one, init_params_kmeans, trapz_weights, normalize_grid_density).
# Every other helper it needs (parameter_names, pack_theta, align_by_location,
# next_power_of_two, aepd_component_scores, make_theory_kde, information_matrix,
# safe_solve) and the scenario definitions are inlined below, so it does NOT
# source scenario1c_asymptotic_normality.R.
#
# Estimator (fully data-driven; the truth is used ONLY to generate data and to
# evaluate errors and coverage, never inside the estimator):
#
#   theta_hat = fixed-order (K0) minimizer of the Hellinger criterion on the
#               active submodel, i.e. the paper's PMHD-MCP-R refit computed with
#               the multi-start fixed_K_fit logic of select_model_stability
#               (lambda = 0, delta_final = max(2/n, 1e-2), k-means starts with a
#               cycled alpha grid, best H2 preferring an all-active solution),
#               evaluated at the inference bandwidth h_n = c * n^{-beta} with
#               beta in (1/(4 abar), 1/2) rather than a Silverman rule.
#
# For each parameter the script reports
#   z_bar          = mean of sqrt(n)(theta_hat - theta*)/SD    (standardized mean)
#   empirical_sd   = sd(theta_hat)                             (Monte-Carlo SD)
#   theoretical_sd = {diag(J_S^{-1})}^{1/2} / sqrt(n)          (plug-in Remark SD)
#   cov95          = empirical coverage of theta_hat +/- z_{0.975} * theoretical_sd
#
# Usage:
#   Rscript scenario1c_asymptotic_normality_final.R --quick
#   Rscript scenario1c_asymptotic_normality_final.R \
#       --scenario=s1d --n=500,1000,2000 --B=200 --cores=12 --bw-constant=2.0
# ============================================================================

script_file <- function() {
  a <- commandArgs(trailingOnly = FALSE); hh <- grep("^--file=", a, value = TRUE)
  if (length(hh)) normalizePath(sub("^--file=", "", hh[1]), mustWork = FALSE)
  else normalizePath("scenario1c_asymptotic_normality_final.R", mustWork = FALSE)
}
script_dir <- dirname(script_file())
core_file <- file.path(script_dir, "simulation_stability_elbow_parallel_main.R")
if (!file.exists(core_file))
  stop("simulation_stability_elbow_parallel_main.R (core algorithm library) not found beside this script.")
suppressPackageStartupMessages(source(core_file))

# ============================================================================
# Inlined helpers (previously in scenario1c_asymptotic_normality.R)
# ============================================================================
parameter_names <- function(K) {
  c(paste0("pi", seq_len(K - 1L)),
    paste0("mu", seq_len(K)),
    paste0("sigma", seq_len(K)),
    paste0("alpha", seq_len(K)),
    paste0("tau", seq_len(K)))
}

pack_theta <- function(par) {
  K <- length(par$pi)
  ans <- c(par$pi[seq_len(K - 1L)], par$mu, par$sigma, par$alpha, par$tau)
  names(ans) <- parameter_names(K)
  ans
}

align_by_location <- function(fit) {
  ord <- order(fit$mu)
  list(pi = fit$pi[ord], mu = fit$mu[ord], sigma = fit$sigma[ord],
       alpha = fit$alpha[ord], tau = fit$tau[ord])
}

next_power_of_two <- function(x) 2^ceiling(log2(max(x, 2)))

# Per-component AEPD score (d log f_k) with respect to (mu, sigma, alpha, tau).
aepd_component_scores <- function(x, mu, sigma, alpha, tau) {
  u <- x - mu
  right <- u >= 0
  q <- ifelse(right, tau, 1 - tau)
  absu <- abs(u)
  base <- q * absu / sigma
  z <- base^alpha
  z_log_base <- numeric(length(base))
  positive <- base > 0
  z_log_base[positive] <- z[positive] * log(base[positive])
  cbind(
    mu = alpha * sign(u) * q^alpha * absu^(alpha - 1) / sigma^alpha,
    sigma = (-1 + alpha * z) / sigma,
    alpha = 1 / alpha + digamma(1 / alpha) / alpha^2 - z_log_base,
    tau = 1 / tau - 1 / (1 - tau) +
      ifelse(right, -alpha * z / tau, alpha * z / (1 - tau))
  )
}

# Epanechnikov-kernel plug-in density on a grid refined so that dx/h is small.
make_theory_kde <- function(x, h, grid_n_min = 2048L, grid_n_max = 131072L,
                            points_per_bandwidth = 10) {
  support_half_width <- sqrt(5) * h
  from <- min(x) - 1.05 * support_half_width
  to <- max(x) + 1.05 * support_half_width
  required <- ceiling((to - from) / (h / points_per_bandwidth)) + 1L
  grid_n <- as.integer(max(grid_n_min, next_power_of_two(required)))
  if (grid_n > grid_n_max) {
    stop("KDE grid requires ", grid_n, " nodes to keep dx/h <= ",
         1 / points_per_bandwidth, ", exceeding the safety limit ", grid_n_max,
         ". Increase kde_grid_n_max; do not silently coarsen the grid.")
  }
  kd <- density(x, bw = h, kernel = "epanechnikov", n = grid_n,
                from = from, to = to, cut = 0)
  wq <- trapz_weights(kd$x)
  grid_step_over_h <- max(diff(kd$x)) / h
  if (!is.finite(grid_step_over_h) ||
      grid_step_over_h > 1 / points_per_bandwidth * 1.001) {
    stop("KDE grid resolution check failed: dx/h = ", signif(grid_step_over_h, 6))
  }
  list(grid = kd$x, wq = wq,
       ghat = normalize_grid_density(kd$y, wq), grid_n = grid_n,
       grid_step_over_h = grid_step_over_h,
       grid_cap_hit = grid_n == grid_n_max)
}

# Hellinger information J_S = 4 integral (nabla sqrt f)(nabla sqrt f)' dx. Under
# correct specification this equals the sandwich covariance of Theorem 3(i).
information_matrix <- function(truth, tail_probability = 1e-10, grid_n = 200001L) {
  K <- length(truth$pi)
  vq <- vapply(seq_len(K), function(k) {
    qgamma(1 - tail_probability, shape = 1 / truth$alpha[k])^(1 / truth$alpha[k])
  }, numeric(1))
  lo <- min(truth$mu - truth$sigma * vq / (1 - truth$tau))
  hi <- max(truth$mu + truth$sigma * vq / truth$tau)
  x <- sort(unique(c(seq(lo, hi, length.out = grid_n), truth$mu)))
  w <- trapz_weights(x)

  comp <- vapply(seq_len(K), function(k) {
    daepd(x, truth$mu[k], truth$sigma[k], truth$alpha[k], truth$tau[k])
  }, numeric(length(x)))
  weighted_comp <- sweep(comp, 2, truth$pi, "*")
  f <- pmax(rowSums(weighted_comp), 1e-300)
  responsibility <- weighted_comp / f

  score <- matrix(0, nrow = length(x), ncol = 5L * K - 1L)
  colnames(score) <- parameter_names(K)
  score[, seq_len(K - 1L)] <- (comp[, seq_len(K - 1L), drop = FALSE] - comp[, K]) / f
  offset <- K - 1L
  for (k in seq_len(K)) {
    cscore <- aepd_component_scores(x, truth$mu[k], truth$sigma[k],
                                    truth$alpha[k], truth$tau[k])
    score[, offset + k]             <- responsibility[, k] * cscore[, "mu"]
    score[, offset + K + k]         <- responsibility[, k] * cscore[, "sigma"]
    score[, offset + 2L * K + k]    <- responsibility[, k] * cscore[, "alpha"]
    score[, offset + 3L * K + k]    <- responsibility[, k] * cscore[, "tau"]
  }
  sw <- sqrt(f * w)
  J <- crossprod(score * sw)
  J <- (J + t(J)) / 2
  attr(J, "captured_mass") <- sum(f * w)
  J
}

safe_solve <- function(A, tolerance = 1e-10) {
  ev <- eigen((A + t(A)) / 2, symmetric = TRUE)
  cutoff <- tolerance * max(ev$values)
  if (min(ev$values) <= cutoff) {
    stop("Information matrix is numerically singular; minimum eigenvalue = ",
         signif(min(ev$values), 4))
  }
  ev$vectors %*% (diag(1 / ev$values, nrow = length(ev$values))) %*% t(ev$vectors)
}

# ============================================================================
# Scenario definitions (three-component AEPD mixtures on mu = (-4, 0, 4))
# ============================================================================
scenario1a <- list(   # Gaussian baseline
  pi = c(0.3, 0.4, 0.3), mu = c(-4, 0, 4),
  sigma = c(0.60, 0.65, 0.60), alpha = c(2.0, 2.0, 2.0), tau = c(0.50, 0.50, 0.50))

scenario1b <- list(   # mild AEPD
  pi = c(0.3, 0.4, 0.3), mu = c(-4, 0, 4),
  sigma = c(0.55, 0.58, 0.55), alpha = c(1.8, 1.7, 1.8), tau = c(0.45, 0.50, 0.55))

scenario1c <- list(   # strong AEPD, original (ill-conditioned) design
  pi = c(0.3, 0.4, 0.3), mu = c(-4, 0, 4),
  sigma = c(0.7, 0.9, 0.7), alpha = c(1.25, 1.15, 1.25), tau = c(0.32, 0.38, 0.44))

scenario1d <- list(   # well-conditioned strong AEPD (paper normality study)
  pi = c(0.3, 0.4, 0.3), mu = c(-4, 0, 4),
  sigma = c(0.5, 0.5, 0.5), alpha = c(1.4, 1.5, 1.3), tau = c(0.45, 0.50, 0.55))

scenario1_settings <- list(s1a = scenario1a, s1b = scenario1b,
                           s1c = scenario1c, s1d = scenario1d)

# ============================================================================
# Estimator and experiment
# ============================================================================
valid_theta <- function(th) {
  all(th$pi > 0.02) && all(th$sigma > 0.2 & th$sigma < 4) &&
    all(th$alpha > 0.5 & th$alpha < 4) &&
    all(th$tau > 0.02 & th$tau < 0.98) && all(is.finite(unlist(th)))
}

# PMHD-MCP-R fixed-order refit (the estimator of Theorem 3(i)); no truth used.
pmhd_mcp_r_refit <- function(x, h, Kc, grid_max, ppb,
                             gamma = 3, refit_nstart = 10L,
                             mm_tol = 1e-5, mm_max_iter = 100L, mm_maxeval = 250L,
                             alpha_lower = 1.0, alpha_upper = 3.0,
                             tau_lower = 0.01, tau_upper = 0.99,
                             active_weight_min = 0.01) {
  kde <- tryCatch(make_theory_kde(x, h, grid_n_max = grid_max, points_per_bandwidth = ppb),
                  error = function(e) NULL)
  if (is.null(kde)) return(NULL)
  n <- length(x)
  delta_final <- max(2 / n, 1e-2)
  alpha_grid  <- c(2.0, 1.5, 1.2, 1.0)
  best_h2 <- Inf; best_fit <- NULL
  best_active_h2 <- Inf; best_active_fit <- NULL
  for (s in seq_len(refit_nstart)) {
    ai <- alpha_grid[((s - 1L) %% length(alpha_grid)) + 1L]
    init_s <- if (s == 1L) NULL
              else tryCatch(init_params_kmeans(x, Kc, alpha_init = ai), error = function(e) NULL)
    f <- tryCatch(
      pmhd_mcp_fit_one(x, Kc, 0, gamma, kde$grid, kde$wq, kde$ghat,
                       tol = mm_tol, max_iter = mm_max_iter,
                       delta_inner = 1e-8, delta_final = delta_final,
                       init = init_s, maxeval = mm_maxeval,
                       alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                       tau_lower = tau_lower, tau_upper = tau_upper),
      error = function(e) NULL)
    if (is.null(f) || is.null(f$Khat) || f$Khat != Kc) next
    if (!all(is.finite(c(f$pi, f$mu, f$sigma, f$alpha, f$tau)))) next
    if (f$H2 < best_h2) { best_h2 <- f$H2; best_fit <- f }
    if (min(f$pi) >= active_weight_min && f$H2 < best_active_h2) {
      best_active_h2 <- f$H2; best_active_fit <- f
    }
  }
  chosen <- if (!is.null(best_active_fit)) best_active_fit else best_fit
  if (is.null(chosen)) return(NULL)
  th <- align_by_location(chosen)
  th <- list(pi = th$pi, mu = th$mu, sigma = th$sigma, alpha = th$alpha, tau = th$tau)
  if (!valid_theta(th)) return(NULL)
  th
}

# One Monte-Carlo replication: fit the estimator and return its root-n error.
process_sample <- function(seed_b, n, h, truth, theta0, Kc,
                           grid_max, ppb, refit_nstart, alpha_lower,
                           mm_tol, mm_max_iter) {
  set.seed(seed_b)
  x <- rmix_aepd(n, truth$pi, truth$mu, truth$sigma, truth$alpha, truth$tau)
  th <- pmhd_mcp_r_refit(x, h, Kc, grid_max, ppb, refit_nstart = refit_nstart,
                         alpha_lower = alpha_lower,
                         mm_tol = mm_tol, mm_max_iter = mm_max_iter)
  if (is.null(th)) return(list(ok = FALSE, reason = "refit"))
  if (max(abs(th$mu - truth$mu)) > 2) return(list(ok = FALSE, reason = "left_branch"))
  list(ok = TRUE, rootn = sqrt(n) * (pack_theta(th) - theta0))
}

# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
run_normality <- function(
    B = 200L, n_values = c(500L, 1000L, 2000L),
    seed = 20240717L, beta = NULL, bandwidth_constant = 2.0, cores = 1L,
    kde_grid_n_max = 262144L, kde_points_per_bandwidth = 10,
    truth_grid_n = 200001L, refit_nstart = 10L, alpha_lower = 1.0,
    mm_tol = 1e-5, mm_max_iter = 100L, scenario = "s1d",
    out_dir = NULL) {

  truth  <- scenario1_settings[[scenario]]
  if (is.null(truth)) stop("Unknown scenario: ", scenario)
  if (is.null(out_dir)) out_dir <- file.path(script_dir, paste0("normality_", scenario))
  theta0 <- pack_theta(truth)
  pnames <- names(theta0)
  Kc     <- length(truth$pi)
  a_bar  <- min(pmin(truth$alpha, 1))
  beta_lower <- 1 / (4 * a_bar)
  if (is.null(beta)) beta <- if (1 / 3 > beta_lower) 1 / 3 else (beta_lower + 0.5) / 2
  if (!(beta > beta_lower && beta < 0.5)) stop("beta must lie in (", signif(beta_lower, 4), ", 0.5).")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat(sprintf("PMHD-MCP-R asymptotic normality  [scenario %s]\n", scenario))
  cat("B =", B, "; n =", paste(n_values, collapse = ", "),
      "; h_n =", signif(bandwidth_constant, 4), "* n^{-", signif(beta, 4), "}\n", sep = "")
  J     <- information_matrix(truth, grid_n = truth_grid_n)
  Sigma <- safe_solve(J)
  se    <- sqrt(diag(Sigma))

  rows <- list(); diag_rows <- list()

  cl <- NULL
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
    # Source THIS self-contained script on each worker: it sources the core
    # library and defines every helper, scenario and estimator function; its
    # main block is guarded by sys.nframe() and does not run on the worker.
    self_file <- script_file()
    parallel::clusterExport(cl, "self_file", envir = environment())
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(source(self_file)))
  }

  for (n in as.integer(n_values)) {
    h <- bandwidth_constant * n^(-beta)
    cat(sprintf("\n n = %d  (h = %.5f)  %d replications ...\n", n, h, B))
    seeds <- seed + n * 1000L + seq_len(B)
    args <- list(n = n, h = h, truth = truth, theta0 = theta0, Kc = Kc,
                 grid_max = kde_grid_n_max, ppb = kde_points_per_bandwidth,
                 refit_nstart = refit_nstart, alpha_lower = alpha_lower,
                 mm_tol = mm_tol, mm_max_iter = mm_max_iter)
    if (is.null(cl)) {
      out <- lapply(seeds, function(sb) do.call(process_sample, c(list(seed_b = sb), args)))
    } else {
      out <- parallel::parLapply(cl, seeds, function(sb) do.call(process_sample, c(list(seed_b = sb), args)))
    }
    ok <- vapply(out, function(z) isTRUE(z$ok), logical(1))
    n_valid <- sum(ok); success_rate <- n_valid / B
    reasons <- table(vapply(out[!ok], function(z) z$reason, character(1)))
    cat(sprintf("   success rate = %.3f (%d / %d valid)\n", success_rate, n_valid, B))
    if (length(reasons)) cat("   failure reasons:", paste(names(reasons), reasons, sep = "=", collapse = ", "), "\n")
    if (n_valid < max(10L, ceiling(0.1 * B))) { cat("   too few valid; skipping n.\n"); next }

    R <- do.call(rbind, lapply(out[ok], `[[`, "rootn"))   # sqrt(n)(theta_hat - theta*)
    Z <- sweep(R, 2, se, "/")                              # standardized

    for (j in seq_along(pnames)) {
      z_bar  <- mean(Z[, j])
      emp_sd <- sd(R[, j]) / sqrt(n)                       # sd(theta_hat)
      th_sd  <- se[j] / sqrt(n)                            # plug-in J^{-1} SD
      cov95  <- mean(abs(Z[, j]) <= qnorm(0.975))
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, parameter = pnames[j], n_valid = n_valid, success_rate = success_rate,
        mean_init = z_bar, empirical_sd = emp_sd, theoretical_sd = th_sd,
        cov95_init = cov95, row.names = NULL)
    }
    diag_rows[[as.character(n)]] <- data.frame(
      n = n, requested = B, valid = n_valid, success_rate = success_rate,
      mean_cov = mean(vapply(seq_along(pnames), function(j) mean(abs(Z[, j]) <= 1.96), numeric(1))))
  }

  tab <- do.call(rbind, rows)
  diagt <- do.call(rbind, diag_rows)
  write.csv(tab,   file.path(out_dir, "normality_coverage.csv"), row.names = FALSE)
  write.csv(diagt, file.path(out_dir, "run_diagnostics.csv"), row.names = FALSE)
  saveRDS(list(config = list(B = B, n_values = n_values, seed = seed, beta = beta,
                             bandwidth_constant = bandwidth_constant, scenario = scenario),
               truth = truth, J = J, Sigma = Sigma, coverage = tab, diagnostics = diagt),
          file.path(out_dir, "normality_results.rds"))

  cat("\nCompleted. Written to:\n ", normalizePath(out_dir), "\n\n")
  print(diagt, row.names = FALSE)
  cat(sprintf("\n%-8s %6s | %7s %10s %10s %6s\n",
              "param", "n", "z_bar", "emp_sd", "theo_sd", "cov95"))
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("%-8s %6d | %7.2f %10.4f %10.4f %6.3f\n",
                tab$parameter[i], tab$n[i], tab$mean_init[i],
                tab$empirical_sd[i], tab$theoretical_sd[i], tab$cov95_init[i]))
  }
  invisible(list(coverage = tab, diagnostics = diagt, J = J, Sigma = Sigma))
}

# ----------------------------------------------------------------------------
parse_cli <- function(args) {
  cfg <- list(B = 200L, n_values = c(500L, 1000L, 2000L), cores = 1L,
              seed = 20240717L, beta = NULL, bandwidth_constant = 2.0,
              refit_nstart = 10L, alpha_lower = 1.0,
              mm_tol = 1e-5, mm_max_iter = 100L, scenario = "s1d",
              out_dir = NULL, quick = FALSE)
  for (a in args) {
    if (a == "--quick") cfg$quick <- TRUE
    else if (grepl("^--B=", a)) cfg$B <- as.integer(sub("^--B=", "", a))
    else if (grepl("^--n=", a)) cfg$n_values <- as.integer(strsplit(sub("^--n=", "", a), ",", fixed = TRUE)[[1]])
    else if (grepl("^--cores=", a)) cfg$cores <- as.integer(sub("^--cores=", "", a))
    else if (grepl("^--seed=", a)) cfg$seed <- as.integer(sub("^--seed=", "", a))
    else if (grepl("^--beta=", a)) cfg$beta <- as.numeric(sub("^--beta=", "", a))
    else if (grepl("^--bw-constant=", a)) cfg$bandwidth_constant <- as.numeric(sub("^--bw-constant=", "", a))
    else if (grepl("^--refit-nstart=", a)) cfg$refit_nstart <- as.integer(sub("^--refit-nstart=", "", a))
    else if (grepl("^--alpha-lower=", a)) cfg$alpha_lower <- as.numeric(sub("^--alpha-lower=", "", a))
    else if (grepl("^--mm-tol=", a)) cfg$mm_tol <- as.numeric(sub("^--mm-tol=", "", a))
    else if (grepl("^--mm-iter=", a)) cfg$mm_max_iter <- as.integer(sub("^--mm-iter=", "", a))
    else if (grepl("^--scenario=", a)) cfg$scenario <- sub("^--scenario=", "", a)
    else if (grepl("^--out=", a)) cfg$out_dir <- sub("^--out=", "", a)
    else stop("Unknown option: ", a)
  }
  if (cfg$quick) { cfg$B <- 20L; cfg$n_values <- 500L; cfg$cores <- 1L }
  cfg
}

if (sys.nframe() == 0L) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  call_args <- list(B = cli$B, n_values = cli$n_values, cores = cli$cores,
                    seed = cli$seed, beta = cli$beta,
                    bandwidth_constant = cli$bandwidth_constant,
                    refit_nstart = cli$refit_nstart, alpha_lower = cli$alpha_lower,
                    mm_tol = cli$mm_tol, mm_max_iter = cli$mm_max_iter,
                    scenario = cli$scenario)
  if (!is.null(cli$out_dir)) call_args$out_dir <- cli$out_dir
  do.call(run_normality, call_args)
}
