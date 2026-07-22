# scenario1c_debiased_normality_implementable.R
# ============================================================================
# Finite-sample check of the Gaussian approximation of Theorem 3(i) for the
# PMHD-MCP-R active-submodel minimizer. This is the experiment behind the
# asymptotic-normality / coverage table in the paper.
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
# The order is fixed at K0, consistent with the conditioning on {Khat = K0} in
# Theorem 3. For each parameter the script reports
#
#   z_bar          = mean of sqrt(n)(theta_hat - theta*)/SD    (standardized mean)
#   empirical_sd   = sd(theta_hat)                             (Monte-Carlo SD)
#   theoretical_sd = {diag(J_S^{-1})}^{1/2} / sqrt(n)          (plug-in Remark SD)
#   cov95          = empirical coverage of theta_hat +/- z_{0.975} * theoretical_sd
#
# No one-step correction, no bootstrap debiasing, no shrinkage: the plug-in
# asymptotic variance J_S^{-1} is validated directly under the well-conditioned
# Scenario 1(d).
#
# Usage:
#   Rscript scenario1c_debiased_normality_implementable.R --quick
#   Rscript scenario1c_debiased_normality_implementable.R \
#       --scenario=s1d --n=500,1000,2000 --B=200 --cores=12 --bw-constant=2.0
# ============================================================================

script_file <- function() {
  a <- commandArgs(trailingOnly = FALSE); hh <- grep("^--file=", a, value = TRUE)
  if (length(hh)) normalizePath(sub("^--file=", "", hh[1]), mustWork = FALSE)
  else normalizePath("scenario1c_debiased_normality_implementable.R", mustWork = FALSE)
}
script_dir <- dirname(script_file())
main_script <- file.path(script_dir, "scenario1c_asymptotic_normality.R")
if (!file.exists(main_script))
  stop("scenario1c_asymptotic_normality.R (function library) not found beside this script.")
suppressPackageStartupMessages(source(main_script))   # main block guarded by sys.nframe()==0

valid_theta <- function(th) {
  all(th$pi > 0.02) && all(th$sigma > 0.2 & th$sigma < 4) &&
    all(th$alpha > 0.5 & th$alpha < 4) &&
    all(th$tau > 0.02 & th$tau < 0.98) && all(is.finite(unlist(th)))
}

# ----------------------------------------------------------------------------
# PMHD-MCP-R fixed-order refit (the estimator of Theorem 3(i)); no truth used.
# Reproduces the multi-start fixed_K_fit logic of select_model_stability at the
# inference bandwidth h.
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# One Monte-Carlo replication: fit the estimator and return its root-n error.
# ----------------------------------------------------------------------------
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
    # Source the full function library on each worker (it in turn sources the
    # core file), so every helper of make_theory_kde / pmhd_mcp_fit_one /
    # rmix_aepd is available; then export the estimator functions defined here.
    main_script_worker <- main_script
    parallel::clusterExport(cl, "main_script_worker", envir = environment())
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(source(main_script_worker)))
    parallel::clusterExport(cl,
      c("pmhd_mcp_r_refit", "valid_theta", "process_sample"),
      envir = globalenv())
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
      z_bar <- mean(Z[, j])
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
