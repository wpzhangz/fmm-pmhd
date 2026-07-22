# Verify the regular-regime asymptotic normality theorem under Scenario 1(c).
#
# Run from a terminal:
#   Rscript scenario1c_asymptotic_normality.R
#   Rscript scenario1c_asymptotic_normality.R --quick
#   Rscript scenario1c_asymptotic_normality.R --B=200 --n=1000,2000,5000 --cores=6
#
# The simulation fixes the correctly selected active order K0=3 and starts
# the optimizer in a neighbourhood of the true parameter.  This targets the
# local active-submodel minimizer whose existence is asserted by Theorem 4,
# rather than mixing its limit law with finite-sample order-selection errors.
# The larger default sample sizes expose the asymptotic trend.  In addition
# to the fitted estimator, the program reports the theorem's influence-
# function linear benchmark and the estimator-minus-linear remainder.

script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  normalizePath("scenario1c_asymptotic_normality.R", mustWork = FALSE)
}

script_dir <- dirname(script_file())
core_file <- file.path(script_dir, "simulation_stability_elbow_parallel_main.R")
if (!file.exists(core_file)) {
  stop("Cannot find simulation_stability_elbow_parallel_main.R beside this script.")
}
suppressPackageStartupMessages(source(core_file))

# Scenario 1(c), copied from paper_scenario1_settings.
scenario1c <- list(
  pi = c(0.3, 0.4, 0.3),
  mu = c(-4, 0, 4),
  sigma = c(0.7, 0.9, 0.7),
  alpha = c(1.25, 1.15, 1.25),
  tau = c(0.32, 0.38, 0.44)
)

# Scenario 1(a): Gaussian baseline. Scales chosen so the mixture peaks form a
# monotone progression 1(a) < 1(b) < 1(d) and all three information matrices
# are well conditioned (cond ~2.6e3, 2.1e3, 1.2e3).
scenario1a <- list(
  pi = c(0.3, 0.4, 0.3),
  mu = c(-4, 0, 4),
  sigma = c(0.60, 0.65, 0.60),
  alpha = c(2.0, 2.0, 2.0),
  tau = c(0.50, 0.50, 0.50)
)

scenario1b <- list(
  pi = c(0.3, 0.4, 0.3),
  mu = c(-4, 0, 4),
  sigma = c(0.55, 0.58, 0.55),
  alpha = c(1.8, 1.7, 1.8),
  tau = c(0.45, 0.50, 0.55)
)

# Scenario 1(d): well-conditioned variant of 1(c) with the SAME location grid
# mu = (-4, 0, 4) and weight structure (0.3, 0.4, 0.3), but narrower components
# (sigma 0.5 -> separation ~8 sigma, so components barely overlap) and shapes
# lifted off the cusp (alpha in [1.3, 1.5]).  This lowers cond(J_S) from ~3400
# to ~1200 and raises n*lambda_min(J_S) at n=500 from 1.1 to 3.8, so the
# asymptotic covariance J_S^{-1} is actually approached at n = 500-2000.
scenario1d <- list(
  pi = c(0.3, 0.4, 0.3),
  mu = c(-4, 0, 4),
  sigma = c(0.5, 0.5, 0.5),
  alpha = c(1.4, 1.5, 1.3),
  tau = c(0.45, 0.50, 0.55)
)

scenario1_settings <- list(s1a = scenario1a, s1b = scenario1b,
                           s1c = scenario1c, s1d = scenario1d)

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

working_from_natural <- function(par) {
  K <- length(par$pi)
  c(log(par$pi[seq_len(K - 1L)] / par$pi[K]),
    par$mu, log(par$sigma), log(par$alpha), qlogis(par$tau))
}

natural_from_working <- function(v, K) {
  i_pi <- seq_len(K - 1L)
  i_mu <- K:(2L * K - 1L)
  i_sigma <- (2L * K):(3L * K - 1L)
  i_alpha <- (3L * K):(4L * K - 1L)
  i_tau <- (4L * K):(5L * K - 1L)
  list(
    pi = softmax(c(v[i_pi], 0)),
    mu = v[i_mu],
    sigma = exp(v[i_sigma]),
    alpha = exp(v[i_alpha]),
    tau = plogis(v[i_tau])
  )
}

softmax <- function(x) {
  y <- exp(x - max(x))
  y / sum(y)
}

make_local_starts <- function(truth, nstart = 3L) {
  starts <- vector("list", nstart)
  starts[[1L]] <- truth
  if (nstart > 1L) {
    for (j in 2:nstart) {
      starts[[j]] <- list(
        pi = softmax(log(truth$pi) + rnorm(length(truth$pi), 0, 0.04)),
        mu = truth$mu + rnorm(length(truth$mu), 0, 0.08),
        sigma = truth$sigma * exp(rnorm(length(truth$sigma), 0, 0.04)),
        alpha = truth$alpha * exp(rnorm(length(truth$alpha), 0, 0.04)),
        tau = plogis(qlogis(truth$tau) + rnorm(length(truth$tau), 0, 0.06))
      )
    }
  }
  starts
}

next_power_of_two <- function(x) 2^ceiling(log2(max(x, 2)))

make_theory_kde <- function(x, h, grid_n_min = 2048L,
                            grid_n_max = 131072L,
                            points_per_bandwidth = 10) {
  # R's Epanechnikov kernel is compactly supported.  The numerical grid is
  # refined as h shrinks so the KDE cusp neighbourhood remains resolved.
  support_half_width <- sqrt(5) * h
  from <- min(x) - 1.05 * support_half_width
  to <- max(x) + 1.05 * support_half_width
  required <- ceiling((to - from) /
                        (h / points_per_bandwidth)) + 1L
  grid_n <- as.integer(max(grid_n_min, next_power_of_two(required)))
  if (grid_n > grid_n_max) {
    stop("KDE grid requires ", grid_n,
         " nodes to keep dx/h <= ", 1 / points_per_bandwidth,
         ", exceeding the safety limit ", grid_n_max,
         ". Increase kde_grid_n_max; do not silently coarsen the grid.")
  }
  kd <- density(x, bw = h, kernel = "epanechnikov", n = grid_n,
                from = from, to = to, cut = 0)
  wq <- trapz_weights(kd$x)
  mass_before_normalizing <- sum(kd$y * wq)
  grid_step_over_h <- max(diff(kd$x)) / h
  if (!is.finite(grid_step_over_h) ||
      grid_step_over_h > 1 / points_per_bandwidth * 1.001) {
    stop("KDE grid resolution check failed: dx/h = ",
         signif(grid_step_over_h, 6))
  }
  list(grid = kd$x, wq = wq,
       ghat = normalize_grid_density(kd$y, wq), grid_n = grid_n,
       grid_n_required = required,
       grid_step_over_h = grid_step_over_h,
       mass_before_normalizing = mass_before_normalizing,
       grid_cap_hit = grid_n == grid_n_max)
}

hellinger_objective_gradient <- function(v, grid, wq, ghat, K) {
  par <- natural_from_working(v, K)
  comp <- vapply(seq_len(K), function(k) {
    daepd(grid, par$mu[k], par$sigma[k], par$alpha[k], par$tau[k])
  }, numeric(length(grid)))
  weighted_comp <- sweep(comp, 2, par$pi, "*")
  f <- pmax(rowSums(weighted_comp), 1e-300)
  responsibility <- weighted_comp / f
  affinity_density <- sqrt(pmax(ghat * f, 0))

  score <- matrix(0, nrow = length(grid), ncol = 5L * K - 1L)
  # Softmax coordinates eta_j=log(pi_j/pi_K).
  score[, seq_len(K - 1L)] <-
    sweep(responsibility[, seq_len(K - 1L), drop = FALSE],
          2, par$pi[seq_len(K - 1L)], "-")
  offset <- K - 1L
  for (k in seq_len(K)) {
    sk <- aepd_component_scores(
      grid, par$mu[k], par$sigma[k], par$alpha[k], par$tau[k]
    )
    score[, offset + k] <- responsibility[, k] * sk[, "mu"]
    score[, offset + K + k] <-
      responsibility[, k] * sk[, "sigma"] * par$sigma[k]
    score[, offset + 2L * K + k] <-
      responsibility[, k] * sk[, "alpha"] * par$alpha[k]
    score[, offset + 3L * K + k] <-
      responsibility[, k] * sk[, "tau"] * par$tau[k] * (1 - par$tau[k])
  }
  aw <- affinity_density * wq
  list(
    objective = 1 - sum(aw),
    gradient = -0.5 * colSums(score * aw)
  )
}

natural_jacobian_from_working <- function(par) {
  K <- length(par$pi)
  d <- 5L * K - 1L
  B <- matrix(0, d, d)
  # The reported natural coordinates contain pi_1,...,pi_{K-1}.
  for (i in seq_len(K - 1L)) {
    for (j in seq_len(K - 1L)) {
      B[i, j] <- par$pi[i] * ((i == j) - par$pi[j])
    }
  }
  i_mu <- K:(2L * K - 1L)
  i_sigma <- (2L * K):(3L * K - 1L)
  i_alpha <- (3L * K):(4L * K - 1L)
  i_tau <- (4L * K):(5L * K - 1L)
  B[i_mu, i_mu] <- diag(K)
  B[i_sigma, i_sigma] <- diag(par$sigma, K)
  B[i_alpha, i_alpha] <- diag(par$alpha, K)
  B[i_tau, i_tau] <- diag(par$tau * (1 - par$tau), K)
  dimnames(B) <- list(parameter_names(K), NULL)
  B
}

finite_difference_gradient_hessian <- function(v, eval_f,
                                                relative_step = 1e-4,
                                                absolute_step = NULL) {
  d <- length(v)
  if (is.null(absolute_step)) {
    absolute_step <- relative_step * pmax(1, abs(v))
  }
  if (length(absolute_step) != d || any(!is.finite(absolute_step)) ||
      any(absolute_step <= 0)) {
    stop("absolute_step must contain one positive finite step per coordinate")
  }
  H <- matrix(NA_real_, d, d)
  for (j in seq_len(d)) {
    eps <- absolute_step[j]
    vp <- vm <- v
    vp[j] <- vp[j] + eps
    vm[j] <- vm[j] - eps
    gp <- eval_f(vp)$gradient
    gm <- eval_f(vm)$gradient
    H[, j] <- (gp - gm) / (2 * eps)
  }
  (H + t(H)) / 2
}

compute_kde_newton_diagnostics <- function(n, truth, grid, wq, ghat, J,
                                           hessian_step = 1e-4) {
  K <- length(truth$pi)
  v0 <- working_from_natural(truth)
  eval_f <- function(v) hellinger_objective_gradient(v, grid, wq, ghat, K)
  at_truth <- eval_f(v0)
  B <- natural_jacobian_from_working(truth)

  # For Q(theta)=1-int sqrt(ghat*f_theta), the population Hessian at the
  # truth is J/4 in natural coordinates. Transform it to working coordinates.
  H_fisher <- crossprod(B, J %*% B) / 4
  fisher_step_v <- -drop(safe_solve(H_fisher) %*% at_truth$gradient)
  fisher_rootn_error <- sqrt(n) * drop(B %*% fisher_step_v)

  coordinate_step <- hessian_step * pmax(1, abs(v0))
  # A location perturbation smaller than the quadrature spacing gives a
  # spurious Hessian at the AEPD cusp. Move each location by at least three
  # grid intervals; the other smooth coordinates retain the relative step.
  location_index <- K:(2L * K - 1L)
  grid_step <- max(diff(grid))
  coordinate_step[location_index] <- pmax(
    coordinate_step[location_index], 3 * grid_step
  )
  H_observed <- finite_difference_gradient_hessian(
    v0, eval_f, relative_step = hessian_step,
    absolute_step = coordinate_step
  )
  eig <- eigen(H_observed, symmetric = TRUE, only.values = TRUE)$values
  observed_step_v <- tryCatch(
    -drop(solve(H_observed, at_truth$gradient)),
    error = function(e) rep(NA_real_, length(v0))
  )
  observed_rootn_error <- sqrt(n) * drop(B %*% observed_step_v)
  names(fisher_rootn_error) <- names(observed_rootn_error) <-
    parameter_names(K)
  list(
    fisher_rootn_error = fisher_rootn_error,
    observed_rootn_error = observed_rootn_error,
    objective_at_truth = at_truth$objective,
    gradient_max_abs_at_truth = max(abs(at_truth$gradient)),
    gradient_l2_at_truth = sqrt(sum(at_truth$gradient^2)),
    observed_hessian_min_eigenvalue = min(eig),
    observed_hessian_max_eigenvalue = max(eig),
    observed_hessian_condition_number =
      max(abs(eig)) / pmax(min(abs(eig)), .Machine$double.eps),
    observed_hessian_positive_definite = min(eig) > 0,
    hessian_location_step_over_grid =
      min(coordinate_step[location_index]) / grid_step
  )
}

direct_refine_active_mhd_wide <- function(init, truth, grid, wq, ghat,
                                     maxeval = 600L,
                                     alpha_lower = 0.55,
                                     alpha_upper = 5,
                                     tau_lower = 0.02,
                                     tau_upper = 0.98) {
  K <- length(truth$pi)
  v0 <- working_from_natural(init)
  # These location boxes isolate the labelled local branch in Scenario 1(c).
  gap <- min(diff(sort(truth$mu)))
  mu_radius <- 0.4 * gap
  lower <- c(rep(-6, K - 1L), truth$mu - mu_radius,
             rep(log(0.08), K), rep(log(alpha_lower), K),
             rep(qlogis(tau_lower), K))
  upper <- c(rep(6, K - 1L), truth$mu + mu_radius,
             rep(log(4), K), rep(log(alpha_upper), K),
             rep(qlogis(tau_upper), K))
  v0 <- pmin(pmax(v0, lower + 1e-8), upper - 1e-8)
  eval_f <- function(v) hellinger_objective_gradient(v, grid, wq, ghat, K)
  opt <- tryCatch(
    nloptr::nloptr(
      x0 = v0, eval_f = eval_f, lb = lower, ub = upper,
      opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-8,
                  ftol_abs = 1e-10, maxeval = maxeval,
                  print_level = 0)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$status < 0 || any(!is.finite(opt$solution))) {
    return(NULL)
  }
  ans <- natural_from_working(opt$solution, K)
  final_eval <- hellinger_objective_gradient(opt$solution, grid, wq, ghat, K)
  ans$H2 <- 2 * max(final_eval$objective, 0)
  ans$optimizer_status <- opt$status
  ans$optimizer_message <- opt$message
  ans$optimizer_iterations <- opt$iterations
  ans$gradient_max_abs <- max(abs(final_eval$gradient))
  ans$boundary_alpha <- any(ans$alpha >= 0.999 * alpha_upper)
  ans$boundary_localization <- any(
    opt$solution - lower <= 1e-5 | upper - opt$solution <= 1e-5
  )
  ans
}

# Oracle local-branch refinement used only to diagnose Theorem 4.  The
# theorem concerns a local minimizer converging to theta*, not an arbitrary
# stationary point of the non-convex mixture objective.  We therefore work
# in root-n, Fisher-whitened coordinates around the truth.  A fixed large
# box in these coordinates contains the Gaussian limit with overwhelming
# probability while excluding remote mixture stationary points.
direct_refine_active_mhd <- function(init, truth, grid, wq, ghat, n, J,
                                     maxeval = 600L,
                                     trust_radius = 8,
                                     hessian_step = 1e-4) {
  K <- length(truth$pi)
  v_truth <- working_from_natural(truth)
  B_truth <- natural_jacobian_from_working(truth)
  B_inv <- solve(B_truth)
  Sigma_v <- B_inv %*% safe_solve(J) %*% t(B_inv)
  Sigma_v <- (Sigma_v + t(Sigma_v)) / 2
  L <- t(chol(Sigma_v))

  eval_v <- function(v) hellinger_objective_gradient(v, grid, wq, ghat, K)
  eval_z <- function(z) {
    v <- v_truth + drop(L %*% z) / sqrt(n)
    ans <- eval_v(v)
    # Multiplying by n leaves the minimizer unchanged and keeps objective
    # differences and gradients on an O(1) scale in root-n coordinates.
    list(
      objective = n * ans$objective,
      gradient = sqrt(n) * drop(crossprod(L, ans$gradient))
    )
  }

  at_truth <- eval_v(v_truth)
  H_fisher_v <- crossprod(B_truth, J %*% B_truth) / 4
  fisher_step_v <- -drop(safe_solve(H_fisher_v) %*% at_truth$gradient)
  fisher_z <- sqrt(n) * drop(solve(L, fisher_step_v))
  if (any(!is.finite(fisher_z))) return(NULL)

  # Backtrack the Fisher step until it is interior and does not increase the
  # numerical Hellinger objective.  This avoids launching L-BFGS across a
  # shallow ridge when the finite-sample Newton step is too long.
  step_scale <- min(1, 0.8 * trust_radius /
                      pmax(max(abs(fisher_z)), .Machine$double.eps))
  z0 <- step_scale * fisher_z
  while (step_scale > 2^-12 &&
         eval_z(z0)$objective > n * at_truth$objective) {
    step_scale <- step_scale / 2
    z0 <- step_scale * fisher_z
  }
  if (eval_z(z0)$objective > n * at_truth$objective) z0[] <- 0

  run_local_opt <- function(start) tryCatch(
    nloptr::nloptr(
      x0 = start, eval_f = eval_z,
      lb = rep(-trust_radius, length(start)),
      ub = rep(trust_radius, length(start)),
      opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-9,
                  ftol_abs = 1e-12, maxeval = maxeval,
                  print_level = 0)
    ),
    error = function(e) NULL
  )
  starts <- list(rep(0, length(z0)), z0 / 2, z0)
  opts <- lapply(starts, run_local_opt)
  valid <- vapply(opts, function(x) {
    # NLopt status -4 means roundoff-limited.  At large n the population
    # target is extremely close to the truth and this benign status is
    # common; retain it provisionally and let the explicit gradient,
    # boundary and curvature checks below decide whether the solution is
    # usable.
    !is.null(x) && (x$status >= 0 || x$status == -4) &&
      all(is.finite(x$solution))
  }, logical(1))
  if (!any(valid)) return(NULL)
  opts <- opts[valid]
  candidate_diagnostics <- lapply(opts, function(x) {
    ez <- eval_z(x$solution)
    c(
      interior = all(abs(x$solution) < 0.999 * trust_radius),
      stationary = max(abs(ez$gradient)) < 1e-4,
      norm = sqrt(sum(x$solution^2)),
      objective = ez$objective
    )
  })
  cd <- do.call(rbind, candidate_diagnostics)
  eligible <- which(cd[, "interior"] > 0.5 & cd[, "stationary"] > 0.5)
  if (!length(eligible)) eligible <- which(cd[, "interior"] > 0.5)
  if (!length(eligible)) eligible <- seq_along(opts)
  # The theorem identifies the local branch approaching theta*.  Among
  # admissible stationary points choose the nearest one, not the globally
  # lowest remote mixture solution.
  chosen <- eligible[which.min(cd[eligible, "norm"])]
  opt <- opts[[chosen]]

  z_hat <- opt$solution
  v_hat <- v_truth + drop(L %*% z_hat) / sqrt(n)
  ans <- natural_from_working(v_hat, K)
  final_eval <- eval_v(v_hat)
  final_eval_z <- eval_z(z_hat)

  # Numerical curvature is a diagnostic, not part of the optimizer.  Move a
  # location by at least three quadrature intervals so that the AEPD cusp is
  # resolved by the finite difference.
  coordinate_step <- hessian_step * pmax(1, abs(v_hat))
  location_index <- K:(2L * K - 1L)
  grid_step <- max(diff(grid))
  coordinate_step[location_index] <- pmax(
    coordinate_step[location_index], 3 * grid_step
  )
  H_hat <- finite_difference_gradient_hessian(
    v_hat, eval_v, relative_step = hessian_step,
    absolute_step = coordinate_step
  )
  eig_hat <- eigen(H_hat, symmetric = TRUE, only.values = TRUE)$values

  ans$H2 <- 2 * max(final_eval$objective, 0)
  ans$optimizer_status <- opt$status
  ans$optimizer_message <- opt$message
  ans$optimizer_iterations <- opt$iterations
  ans$gradient_max_abs <- max(abs(final_eval$gradient))
  ans$local_gradient_max_abs <- max(abs(final_eval_z$gradient))
  ans$boundary_alpha <- FALSE
  ans$boundary_localization <- any(abs(z_hat) >= 0.999 * trust_radius)
  ans$trust_coordinate_max_abs <- max(abs(z_hat))
  ans$trust_coordinate_l2 <- sqrt(sum(z_hat^2))
  ans$fisher_start_max_abs <- max(abs(fisher_z))
  ans$fisher_start_scale <- step_scale
  ans$local_starts_attempted <- length(starts)
  ans$local_candidates_valid <- length(opts)
  ans$local_candidate_selected <- chosen
  ans$distance_from_fisher_l2 <- sqrt(sum((z_hat - fisher_z)^2))
  ans$objective_at_truth <- at_truth$objective
  ans$objective_decrease_from_truth <- at_truth$objective - final_eval$objective
  ans$hessian_min_eigenvalue <- min(eig_hat)
  ans$hessian_max_eigenvalue <- max(eig_hat)
  ans$hessian_positive_definite <- min(eig_hat) > 0
  ans
}

fit_active_mhd <- function(x, truth, h, lambda_n, gamma = 3,
                           nstart = 3L, max_iter = 150L,
                           maxeval = 600L, direct_refine = TRUE,
                           J = NULL, local_trust_radius = 8,
                           newton_hessian_step = 1e-4,
                           kde_grid_n_max = 131072L,
                           kde_points_per_bandwidth = 10,
                           kde = NULL) {
  if (is.null(kde)) {
    kde <- make_theory_kde(
      x, h, grid_n_max = kde_grid_n_max,
      points_per_bandwidth = kde_points_per_bandwidth
    )
  }
  # In the oracle normality experiment the active order and truth-centred
  # local chart are fixed by design.  The MCP-MM stage is irrelevant on the
  # flat tail and can send the calculation to another mixture branch before
  # the local refinement starts.  Bypass it here and optimize the active
  # Hellinger objective directly in the theorem's local chart.
  if (isTRUE(direct_refine) && !is.null(J)) {
    ans <- direct_refine_active_mhd(
      truth, truth, kde$grid, kde$wq, kde$ghat,
      n = length(x), J = J, maxeval = maxeval,
      trust_radius = local_trust_radius,
      hessian_step = newton_hessian_step
    )
    if (is.null(ans)) return(NULL)
    ans$grid_n <- kde$grid_n
    ans$grid_n_required <- kde$grid_n_required
    ans$grid_step_over_h <- kde$grid_step_over_h
    ans$kde_mass_before_normalizing <- kde$mass_before_normalizing
    ans$grid_cap_hit <- kde$grid_cap_hit
    ans$flat_tail <- all(ans$pi > gamma * lambda_n)
    return(ans)
  }
  starts <- make_local_starts(truth, nstart)
  fits <- lapply(starts, function(init) {
    tryCatch(
      pmhd_mcp_fit_one(
        x = x, K_init = length(truth$pi), lambda = lambda_n, gamma = gamma,
        grid = kde$grid, wq = kde$wq, ghat = kde$ghat,
        tol = 1e-7, max_iter = max_iter, delta_inner = 1e-10,
        delta_final = 0, init = init, maxeval = maxeval,
        alpha_lower = 0.55, alpha_upper = 3,
        tau_lower = 0.02, tau_upper = 0.98
      ),
      error = function(e) NULL
    )
  })
  valid <- vapply(fits, function(z) {
    !is.null(z) && z$Khat == length(truth$pi) && is.finite(z$H2) &&
      all(is.finite(c(z$pi, z$mu, z$sigma, z$alpha, z$tau)))
  }, logical(1))
  if (!any(valid)) return(NULL)
  fits <- fits[valid]
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "H2"))]]
  ans <- align_by_location(best)
  if (isTRUE(direct_refine)) {
    if (is.null(J)) stop("J is required for local-branch refinement")
    refined <- direct_refine_active_mhd(
      ans, truth, kde$grid, kde$wq, kde$ghat,
      n = length(x), J = J, maxeval = maxeval,
      trust_radius = local_trust_radius,
      hessian_step = newton_hessian_step
    )
    if (is.null(refined)) return(NULL)
    ans <- refined
  } else {
    ans$H2 <- best$H2
    ans$optimizer_status <- NA_integer_
    ans$optimizer_message <- "MM only"
    ans$optimizer_iterations <- NA_integer_
    ans$gradient_max_abs <- NA_real_
    ans$boundary_alpha <- any(ans$alpha >= 0.999 * 3)
    ans$boundary_localization <- NA
  }
  ans$grid_n <- kde$grid_n
  ans$grid_n_required <- kde$grid_n_required
  ans$grid_step_over_h <- kde$grid_step_over_h
  ans$kde_mass_before_normalizing <- kde$mass_before_normalizing
  ans$grid_cap_hit <- kde$grid_cap_hit
  ans$flat_tail <- all(ans$pi > gamma * lambda_n)
  ans
}

# One Fisher-scoring update of the KDE score equation, using the converged
# PMHD estimate only as a stable local-branch pilot.  The update is computed
# in the natural simplex chart.  This avoids exponentiating a non-negligible
# working-coordinate step, which was unstable in weakly identified shape
# directions.  A spectral floor, a Fisher-norm trust bound and step halving
# keep the update on the selected local branch.
#
# The score moment uses the KDE, whereas the information is evaluated under
# the fitted parametric model.  This gives the usual one-step cancellation
#   theta_os - theta0 = J(theta0)^{-1} int ghat score(theta0) + o_p(n^-1/2)
# whenever the pilot is sufficiently accurate.
kde_score_one_step <- function(pilot, grid, wq, ghat, n,
                               ridge_relative = 1e-3,
                               max_rootn_fisher_step = 6,
                               max_step_halving = 30L) {
  K <- length(pilot$pi)
  theta_pilot <- pack_theta(pilot)
  score_natural <- mixture_score_matrix(grid, pilot)

  comp <- vapply(seq_len(K), function(k) {
    daepd(grid, pilot$mu[k], pilot$sigma[k], pilot$alpha[k], pilot$tau[k])
  }, numeric(length(grid)))
  f_pilot <- pmax(drop(comp %*% pilot$pi), 1e-300)
  # The KDE grid is data-dependent and may omit negligible parametric tails.
  # Do not renormalize f_pilot: the captured-mass diagnostic records this
  # numerical truncation explicitly.
  # Center by the model score on the same finite grid.  Its full-support
  # integral is exactly zero, so this is a numerical control variate that
  # removes grid-tail truncation without changing the first-order equation.
  moment_natural <- colSums(
    score_natural * ((ghat - f_pilot) * wq)
  )
  information_natural <- crossprod(
    score_natural * sqrt(f_pilot * wq)
  )
  information_natural <-
    (information_natural + t(information_natural)) / 2
  ev <- eigen(information_natural, symmetric = TRUE)
  cutoff <- ridge_relative * max(ev$values)
  if (!is.finite(cutoff) || cutoff <= 0) stop("invalid one-step information")
  stabilized_values <- pmax(ev$values, cutoff)
  inverse_information <- ev$vectors %*%
    (diag(1 / stabilized_values, nrow = length(stabilized_values))) %*%
    t(ev$vectors)
  delta <- drop(inverse_information %*% moment_natural)
  fisher_norm <- sqrt(max(0, n * drop(
    crossprod(delta, information_natural %*% delta)
  )))
  clipped <- FALSE
  if (is.finite(max_rootn_fisher_step) &&
      is.finite(fisher_norm) && fisher_norm > max_rootn_fisher_step) {
    delta <- delta * max_rootn_fisher_step / fisher_norm
    fisher_norm <- max_rootn_fisher_step
    clipped <- TRUE
  }
  feasible <- function(theta) {
    par <- unpack_theta_natural(theta, K)
    all(is.finite(theta)) && all(par$pi > 1e-6) &&
      all(par$sigma > 0.05 & par$sigma < 10) &&
      all(par$alpha > 0.55 & par$alpha < 5) &&
      all(par$tau > 0.02 & par$tau < 0.98) &&
      max(abs(par$mu - pilot$mu)) <= 2
  }
  step_scale <- 1
  theta_updated <- theta_pilot + delta
  halvings <- 0L
  while (!feasible(theta_updated) && halvings < max_step_halving) {
    step_scale <- step_scale / 2
    theta_updated <- theta_pilot + step_scale * delta
    halvings <- halvings + 1L
  }
  if (!feasible(theta_updated)) stop("one-step update left feasible local chart")
  updated <- unpack_theta_natural(theta_updated, K)
  list(
    estimate = updated,
    theta = theta_updated,
    moment_max_abs = max(abs(moment_natural)),
    moment_l2 = sqrt(sum(moment_natural^2)),
    rootn_fisher_step = fisher_norm * step_scale,
    information_min_eigenvalue = min(ev$values),
    information_condition_number = max(ev$values) / min(ev$values),
    information_ridge_used = any(ev$values < cutoff),
    model_mass_on_kde_grid = sum(f_pilot * wq),
    step_clipped = clipped,
    step_scale = step_scale,
    step_halvings = halvings
  )
}

unpack_theta_natural <- function(v, K) {
  ip <- seq_len(K - 1L)
  ans <- list(
    pi = c(v[ip], 1 - sum(v[ip])),
    mu = v[K:(2L * K - 1L)],
    sigma = v[(2L * K):(3L * K - 1L)],
    alpha = v[(3L * K):(4L * K - 1L)],
    tau = v[(4L * K):(5L * K - 1L)]
  )
  ans
}

# Feasible per-sample parametric-bootstrap centering of the score one-step.
# No truth is used.  The optional shrinkage factor estimates the signal-to-
# Monte-Carlo-noise ratio of the bootstrap centering term; unlike the earlier
# diagnostic script it is not tuned against the outer simulation error.
kde_score_bootstrap_centering <- function(
    estimate, n, h, Bboot, grid_n_max, points_per_bandwidth,
    ridge_relative, max_rootn_fisher_step,
    min_valid_fraction = 0.7) {
  K <- length(estimate$pi)
  theta_center <- pack_theta(estimate)
  draws <- matrix(NA_real_, Bboot, length(theta_center))
  for (j in seq_len(Bboot)) {
    xb <- rmix_aepd(n, estimate$pi, estimate$mu, estimate$sigma,
                    estimate$alpha, estimate$tau)
    kde_b <- tryCatch(make_theory_kde(
      xb, h, grid_n_max = grid_n_max,
      points_per_bandwidth = points_per_bandwidth
    ), error = function(e) NULL)
    if (is.null(kde_b)) next
    os_b <- tryCatch(kde_score_one_step(
      estimate, kde_b$grid, kde_b$wq, kde_b$ghat, n,
      ridge_relative = ridge_relative,
      max_rootn_fisher_step = max_rootn_fisher_step
    ), error = function(e) NULL)
    if (!is.null(os_b)) draws[j, ] <- sqrt(n) * (os_b$theta - theta_center)
  }
  keep <- complete.cases(draws)
  need <- max(2L, ceiling(min_valid_fraction * Bboot))
  if (sum(keep) < need) stop("too few valid score bootstrap draws")
  draws <- draws[keep, , drop = FALSE]
  bias <- colMeans(draws)
  bias_mc_variance <- apply(draws, 2, var) / nrow(draws)
  shrinkage <- bias^2 / (bias^2 + bias_mc_variance)
  shrinkage[!is.finite(shrinkage)] <- 0
  shrinkage <- pmin(1, pmax(0, shrinkage))

  apply_correction <- function(multiplier) {
    direction <- multiplier * bias / sqrt(n)
    scale <- 1
    feasible <- function(theta) {
      par <- unpack_theta_natural(theta, K)
      all(is.finite(theta)) && all(par$pi > 1e-6) &&
        all(par$sigma > 0.05 & par$sigma < 10) &&
        all(par$alpha > 0.55 & par$alpha < 5) &&
        all(par$tau > 0.02 & par$tau < 0.98)
    }
    candidate <- theta_center - direction
    while (!feasible(candidate) && scale > 2^-30) {
      scale <- scale / 2
      candidate <- theta_center - scale * direction
    }
    if (!feasible(candidate)) stop("bootstrap correction left parameter space")
    list(theta = candidate, scale = scale)
  }
  full <- apply_correction(rep(1, length(bias)))
  shrunk <- apply_correction(shrinkage)
  list(
    bias_rootn = bias,
    bias_mcse_rootn = sqrt(bias_mc_variance),
    shrinkage = shrinkage,
    full_theta = full$theta,
    shrunk_theta = shrunk$theta,
    full_scale = full$scale,
    shrunk_scale = shrunk$scale,
    valid = nrow(draws)
  )
}

simulate_one_normality <- function(b, n, truth, beta, bandwidth_constant,
                                   lambda_constant, lambda_power, gamma,
                                   nstart, max_iter, maxeval, Sigma, J,
                                   kde_grid_n_max,
                                   kde_points_per_bandwidth,
                                   newton_hessian_step,
                                   local_trust_radius,
                                   score_ridge_relative,
                                   score_max_rootn_step,
                                   score_bootstrap_B,
                                   score_bootstrap_min_valid_fraction) {
  x <- rmix_aepd(n, truth$pi, truth$mu, truth$sigma,
                 truth$alpha, truth$tau)
  score_at_truth <- mixture_score_matrix(x, truth)
  linear_rootn_error <- drop(Sigma %*% (colSums(score_at_truth) / sqrt(n)))
  names(linear_rootn_error) <- colnames(score_at_truth)
  na_vector <- setNames(rep(NA_real_, length(linear_rootn_error)),
                        names(linear_rootn_error))
  empty_newton <- list(
    fisher_rootn_error = na_vector,
    observed_rootn_error = na_vector,
    objective_at_truth = NA_real_, gradient_max_abs_at_truth = NA_real_,
    gradient_l2_at_truth = NA_real_,
    observed_hessian_min_eigenvalue = NA_real_,
    observed_hessian_max_eigenvalue = NA_real_,
    observed_hessian_condition_number = NA_real_,
    observed_hessian_positive_definite = FALSE,
    hessian_location_step_over_grid = NA_real_
  )
  h <- bandwidth_constant * n^(-beta)
  lambda_n <- lambda_constant * n^(-lambda_power)
  kde <- tryCatch(
    make_theory_kde(
      x, h, grid_n_max = kde_grid_n_max,
      points_per_bandwidth = kde_points_per_bandwidth
    ),
    error = function(e) NULL
  )
  if (is.null(kde)) {
    return(c(list(ok = FALSE, b = b, message = "KDE construction failed",
                  linear_rootn_error = linear_rootn_error,
                  kde_linear_rootn_error = na_vector), empty_newton))
  }
  # Exact finite-sample version of the leading KDE term in the proof.  Under
  # correct specification u(x)=nabla sqrt(p0)/sqrt(p0)=score(x)/2 and
  # integral p0(x) score(x) dx=0, hence
  #   2 J^{-1} barZ_n = J^{-1} sqrt(n) integral ghat(x) score(x) dx.
  score_on_grid <- mixture_score_matrix(kde$grid, truth)
  kde_score_mean <- colSums(
    score_on_grid * (kde$ghat * kde$wq)
  )
  kde_linear_rootn_error <-
    sqrt(n) * drop(Sigma %*% kde_score_mean)
  names(kde_linear_rootn_error) <- colnames(score_on_grid)
  newton <- tryCatch(
    compute_kde_newton_diagnostics(
      n, truth, kde$grid, kde$wq, kde$ghat, J,
      hessian_step = newton_hessian_step
    ),
    error = function(e) NULL
  )
  if (is.null(newton)) newton <- empty_newton
  fail <- function(message, fit = NULL) {
    ans <- c(list(ok = FALSE, b = b, message = message,
                  linear_rootn_error = linear_rootn_error,
                  kde_linear_rootn_error = kde_linear_rootn_error), newton)
    if (!is.null(fit)) {
      ans$candidate_theta <- pack_theta(fit)
      ans$candidate_gradient_max_abs <- fit$gradient_max_abs
      ans$candidate_local_gradient_max_abs <- fit$local_gradient_max_abs
      ans$candidate_trust_coordinate_max_abs <-
        fit$trust_coordinate_max_abs
      ans$candidate_hessian_min_eigenvalue <-
        fit$hessian_min_eigenvalue
      ans$candidate_objective_decrease_from_truth <-
        fit$objective_decrease_from_truth
    }
    ans
  }
  fit <- fit_active_mhd(
    x, truth, h, lambda_n, gamma, nstart, max_iter, maxeval,
    kde_grid_n_max = kde_grid_n_max,
    kde_points_per_bandwidth = kde_points_per_bandwidth,
    kde = kde, J = J, local_trust_radius = local_trust_radius,
    newton_hessian_step = newton_hessian_step
  )
  if (is.null(fit)) {
    return(fail("fit failed"))
  }
  # The separated Scenario 1(c) locations make sorting an unambiguous label map.
  if (max(abs(fit$mu - truth$mu)) > 2) {
    return(fail("left local branch", fit))
  }
  if (!isTRUE(fit$flat_tail)) {
    return(fail("left MCP flat-tail neighbourhood", fit))
  }
  if (isTRUE(fit$boundary_localization)) {
    return(fail("refinement hit a localization/parameter boundary", fit))
  }
  if (!is.finite(fit$gradient_max_abs) || fit$gradient_max_abs > 5e-3) {
    return(fail("joint objective did not reach gradient tolerance", fit))
  }
  if (!is.finite(fit$local_gradient_max_abs) ||
      fit$local_gradient_max_abs > 1e-4) {
    return(fail("local root-n objective did not reach gradient tolerance",
                fit))
  }
  if (!isTRUE(fit$hessian_positive_definite)) {
    return(fail("local objective Hessian is not positive definite", fit))
  }
  if (!is.finite(fit$objective_decrease_from_truth) ||
      fit$objective_decrease_from_truth < -1e-10) {
    return(fail("local objective did not improve on the truth", fit))
  }
  score_one_step <- tryCatch(
    kde_score_one_step(
      fit, kde$grid, kde$wq, kde$ghat, n,
      ridge_relative = score_ridge_relative,
      max_rootn_fisher_step = score_max_rootn_step
    ),
    error = function(e) NULL
  )
  if (is.null(score_one_step)) {
    score_one_step <- list(
      theta = setNames(rep(NA_real_, length(pack_theta(truth))),
                       names(pack_theta(truth))),
      moment_max_abs = NA_real_, moment_l2 = NA_real_,
      rootn_fisher_step = NA_real_,
      information_min_eigenvalue = NA_real_,
      information_condition_number = NA_real_,
      information_ridge_used = NA,
      model_mass_on_kde_grid = NA_real_, step_clipped = NA,
      step_scale = NA_real_, step_halvings = NA_integer_
    )
  }
  score_bootstrap <- NULL
  if (score_bootstrap_B > 0L && all(is.finite(score_one_step$theta))) {
    score_bootstrap <- tryCatch(kde_score_bootstrap_centering(
      unpack_theta_natural(score_one_step$theta, length(truth$pi)),
      n = n, h = h, Bboot = score_bootstrap_B,
      grid_n_max = kde_grid_n_max,
      points_per_bandwidth = kde_points_per_bandwidth,
      ridge_relative = score_ridge_relative,
      max_rootn_fisher_step = score_max_rootn_step,
      min_valid_fraction = score_bootstrap_min_valid_fraction
    ), error = function(e) NULL)
  }
  if (is.null(score_bootstrap)) {
    score_bootstrap <- list(
      bias_rootn = na_vector, bias_mcse_rootn = na_vector,
      shrinkage = na_vector, full_theta = na_vector,
      shrunk_theta = na_vector, full_scale = NA_real_,
      shrunk_scale = NA_real_, valid = 0L
    )
  }
  c(list(ok = TRUE, b = b, theta = pack_theta(fit), H2 = fit$H2,
       kde_score_one_step_theta = score_one_step$theta,
       kde_score_bootstrap_full_theta = score_bootstrap$full_theta,
       kde_score_bootstrap_shrunk_theta = score_bootstrap$shrunk_theta,
       kde_score_bootstrap_bias_rootn = score_bootstrap$bias_rootn,
       kde_score_bootstrap_bias_mcse_rootn =
         score_bootstrap$bias_mcse_rootn,
       kde_score_bootstrap_shrinkage = score_bootstrap$shrinkage,
       kde_score_bootstrap_valid = score_bootstrap$valid,
       kde_score_bootstrap_full_scale = score_bootstrap$full_scale,
       kde_score_bootstrap_shrunk_scale = score_bootstrap$shrunk_scale,
       kde_score_moment_max_abs = score_one_step$moment_max_abs,
       kde_score_moment_l2 = score_one_step$moment_l2,
       kde_score_rootn_fisher_step = score_one_step$rootn_fisher_step,
       kde_score_information_min_eigenvalue =
         score_one_step$information_min_eigenvalue,
       kde_score_information_condition_number =
         score_one_step$information_condition_number,
       kde_score_information_ridge_used =
         score_one_step$information_ridge_used,
       kde_score_model_mass_on_grid = score_one_step$model_mass_on_kde_grid,
       kde_score_step_clipped = score_one_step$step_clipped,
       kde_score_step_scale = score_one_step$step_scale,
       kde_score_step_halvings = score_one_step$step_halvings,
       flat_tail = fit$flat_tail, grid_n = fit$grid_n,
       grid_n_required = fit$grid_n_required,
       grid_step_over_h = fit$grid_step_over_h,
       kde_mass_before_normalizing = fit$kde_mass_before_normalizing,
       grid_cap_hit = fit$grid_cap_hit,
       h = h, lambda = lambda_n,
       optimizer_status = fit$optimizer_status,
       optimizer_iterations = fit$optimizer_iterations,
       gradient_max_abs = fit$gradient_max_abs,
       local_gradient_max_abs = fit$local_gradient_max_abs,
       trust_coordinate_max_abs = fit$trust_coordinate_max_abs,
       trust_coordinate_l2 = fit$trust_coordinate_l2,
       fisher_start_max_abs = fit$fisher_start_max_abs,
       fisher_start_scale = fit$fisher_start_scale,
       local_starts_attempted = fit$local_starts_attempted,
       local_candidates_valid = fit$local_candidates_valid,
       local_candidate_selected = fit$local_candidate_selected,
       distance_from_fisher_l2 = fit$distance_from_fisher_l2,
       objective_decrease_from_truth = fit$objective_decrease_from_truth,
       hessian_min_eigenvalue = fit$hessian_min_eigenvalue,
       hessian_positive_definite = fit$hessian_positive_definite,
       boundary_alpha = fit$boundary_alpha,
       boundary_localization = fit$boundary_localization,
       linear_rootn_error = linear_rootn_error,
       kde_linear_rootn_error = kde_linear_rootn_error), newton)
}

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

mixture_score_matrix <- function(x, par) {
  K <- length(par$pi)
  comp <- vapply(seq_len(K), function(k) {
    daepd(x, par$mu[k], par$sigma[k], par$alpha[k], par$tau[k])
  }, numeric(length(x)))
  weighted_comp <- sweep(comp, 2, par$pi, "*")
  f <- pmax(rowSums(weighted_comp), 1e-300)
  responsibility <- weighted_comp / f
  score <- matrix(0, nrow = length(x), ncol = 5L * K - 1L)
  colnames(score) <- parameter_names(K)
  score[, seq_len(K - 1L)] <-
    (comp[, seq_len(K - 1L), drop = FALSE] - comp[, K]) / f
  offset <- K - 1L
  for (k in seq_len(K)) {
    sk <- aepd_component_scores(
      x, par$mu[k], par$sigma[k], par$alpha[k], par$tau[k]
    )
    score[, offset + k] <- responsibility[, k] * sk[, "mu"]
    score[, offset + K + k] <- responsibility[, k] * sk[, "sigma"]
    score[, offset + 2L * K + k] <- responsibility[, k] * sk[, "alpha"]
    score[, offset + 3L * K + k] <- responsibility[, k] * sk[, "tau"]
  }
  score
}

information_matrix <- function(truth, tail_probability = 1e-10,
                               grid_n = 200001L) {
  K <- length(truth$pi)
  # Tail limits are obtained from the exact gamma representation used by raepd.
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
  score[, seq_len(K - 1L)] <-
    (comp[, seq_len(K - 1L), drop = FALSE] - comp[, K]) / f

  offset <- K - 1L
  component_score <- vector("list", K)
  for (k in seq_len(K)) {
    component_score[[k]] <- aepd_component_scores(
      x, truth$mu[k], truth$sigma[k], truth$alpha[k], truth$tau[k])
    score[, offset + k] <- responsibility[, k] * component_score[[k]][, "mu"]
    score[, offset + K + k] <- responsibility[, k] * component_score[[k]][, "sigma"]
    score[, offset + 2L * K + k] <- responsibility[, k] * component_score[[k]][, "alpha"]
    score[, offset + 3L * K + k] <- responsibility[, k] * component_score[[k]][, "tau"]
  }

  # Under correct specification this is both the Fisher information and
  # J_S = 4 integral (nabla sqrt(f))(nabla sqrt(f))' dx in Theorem 4.
  sw <- sqrt(f * w)
  J <- crossprod(score * sw)
  J <- (J + t(J)) / 2
  attr(J, "captured_mass") <- sum(f * w)
  attr(J, "limits") <- c(lo, hi)
  J
}

safe_solve <- function(A, tolerance = 1e-10) {
  ev <- eigen((A + t(A)) / 2, symmetric = TRUE)
  cutoff <- tolerance * max(ev$values)
  if (min(ev$values) <= cutoff) {
    stop("Information matrix is numerically singular; minimum eigenvalue = ",
         signif(min(ev$values), 4))
  }
  ev$vectors %*% (diag(1 / ev$values, nrow = length(ev$values))) %*%
    t(ev$vectors)
}

true_mixture_density <- function(x, truth) {
  Reduce("+", lapply(seq_along(truth$pi), function(k) {
    truth$pi[k] * daepd(
      x, truth$mu[k], truth$sigma[k], truth$alpha[k], truth$tau[k]
    )
  }))
}

make_population_smoothed_density <- function(
    truth, h, tail_probability = 1e-10,
    points_per_bandwidth = 10, kernel_nodes = 401L,
    grid_n_max = 131072L) {
  K <- length(truth$pi)
  kernel_nodes <- as.integer(kernel_nodes)
  if (kernel_nodes < 41L) stop("population kernel_nodes must be at least 41")
  if (kernel_nodes %% 2L == 0L) kernel_nodes <- kernel_nodes + 1L

  vq <- vapply(seq_len(K), function(k) {
    qgamma(1 - tail_probability, shape = 1 / truth$alpha[k])^(
      1 / truth$alpha[k]
    )
  }, numeric(1))
  lo <- min(truth$mu - truth$sigma * vq / (1 - truth$tau)) - sqrt(5) * h
  hi <- max(truth$mu + truth$sigma * vq / truth$tau) + sqrt(5) * h
  required <- ceiling((hi - lo) / (h / points_per_bandwidth)) + 1L
  grid_n <- as.integer(next_power_of_two(required))
  if (grid_n > grid_n_max) {
    stop("Population-smoothed grid requires ", grid_n,
         " nodes, exceeding grid_n_max=", grid_n_max)
  }
  grid <- seq(lo, hi, length.out = grid_n)
  wq <- trapz_weights(grid)

  u <- seq(-sqrt(5), sqrt(5), length.out = kernel_nodes)
  wu <- trapz_weights(u)
  kernel_value <- 3 / (4 * sqrt(5)) * pmax(1 - u^2 / 5, 0)
  kernel_weight <- kernel_value * wu
  kernel_weight <- kernel_weight / sum(kernel_weight)
  p0h <- numeric(length(grid))
  for (j in seq_along(u)) {
    p0h <- p0h + kernel_weight[j] *
      true_mixture_density(grid - h * u[j], truth)
  }
  mass_before_normalizing <- sum(p0h * wq)
  p0h <- normalize_grid_density(p0h, wq)
  list(
    grid = grid, wq = wq, p0h = p0h,
    grid_n = grid_n, grid_n_required = required,
    grid_step_over_h = max(diff(grid)) / h,
    mass_before_normalizing = mass_before_normalizing,
    kernel_nodes = kernel_nodes
  )
}

compute_population_smoothed_target <- function(
    n, truth, beta, bandwidth_constant, Sigma,
    maxeval = 2000L, kernel_nodes = 401L,
    grid_n_max = 131072L, points_per_bandwidth = 10,
    local_trust_radius = 8) {
  h <- bandwidth_constant * n^(-beta)
  smoothed <- make_population_smoothed_density(
    truth, h, points_per_bandwidth = points_per_bandwidth,
    kernel_nodes = kernel_nodes, grid_n_max = grid_n_max
  )
  # This is a deterministic population target, not a random root-n local
  # estimator.  Use the ordinary truth-initialized refinement; imposing the
  # sample-size-dependent Fisher trust chart can become numerically
  # degenerate when the smoothing bias is extremely small.
  fit_h <- direct_refine_active_mhd_wide(
    init = truth, truth = truth,
    grid = smoothed$grid, wq = smoothed$wq, ghat = smoothed$p0h,
    maxeval = maxeval
  )
  if (is.null(fit_h) || isTRUE(fit_h$boundary_localization) ||
      !is.finite(fit_h$gradient_max_abs) ||
      fit_h$gradient_max_abs > 1e-4) {
    stop("Population-smoothed optimization failed at n=", n)
  }
  theta_true <- pack_theta(truth)
  theta_h <- pack_theta(fit_h)
  bias <- theta_h - theta_true
  data.frame(
    n = n, bandwidth = h, parameter = names(theta_true),
    theta_true = as.numeric(theta_true), theta_h = as.numeric(theta_h),
    smoothing_bias = as.numeric(bias),
    root_n_smoothing_bias = sqrt(n) * as.numeric(bias),
    standardized_smoothing_bias =
      sqrt(n) * as.numeric(bias) / sqrt(diag(Sigma)),
    population_H2 = fit_h$H2,
    optimizer_status = fit_h$optimizer_status,
    gradient_max_abs = fit_h$gradient_max_abs,
    boundary_localization = fit_h$boundary_localization,
    grid_n = smoothed$grid_n,
    grid_n_required = smoothed$grid_n_required,
    grid_step_over_h = smoothed$grid_step_over_h,
    mass_before_normalizing = smoothed$mass_before_normalizing,
    kernel_nodes = smoothed$kernel_nodes,
    stringsAsFactors = FALSE
  )
}

normal_summary <- function(z, n, names_z) {
  B <- nrow(z)
  covered_95 <- colMeans(abs(z) <= qnorm(0.975))
  data.frame(
    n = n,
    parameter = names_z,
    replications = B,
    mean = colMeans(z),
    sd = apply(z, 2, sd),
    skewness = apply(z, 2, function(v) mean((v - mean(v))^3) / sd(v)^3),
    excess_kurtosis = apply(z, 2, function(v) mean((v - mean(v))^4) / sd(v)^4 - 3),
    qq_correlation = apply(z, 2, function(v) {
      cor(sort(v), qnorm(ppoints(length(v))))
    }),
    coverage_95 = covered_95,
    coverage_95_mcse = sqrt(covered_95 * (1 - covered_95) / B),
    stringsAsFactors = FALSE
  )
}

add_standard_deviation_comparison <- function(summary_table, Sigma,
                                              parameter_order) {
  theoretical_root_n_sd <- setNames(
    sqrt(diag(Sigma)), parameter_order
  )
  scale_root_n <- theoretical_root_n_sd[summary_table$parameter]
  summary_table$standardized_sd <- summary_table$sd
  summary_table$empirical_root_n_sd <- summary_table$sd * scale_root_n
  summary_table$theoretical_root_n_sd <- scale_root_n
  summary_table$empirical_sd <-
    summary_table$empirical_root_n_sd / sqrt(summary_table$n)
  summary_table$theoretical_sd <-
    summary_table$theoretical_root_n_sd / sqrt(summary_table$n)
  summary_table$empirical_to_theoretical_sd_ratio <-
    summary_table$empirical_sd / summary_table$theoretical_sd
  summary_table
}

plot_marginal_qq <- function(z_by_n, file) {
  grDevices::pdf(file, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  for (nm in names(z_by_n)) {
    z <- z_by_n[[nm]]
    par(mfrow = c(4, 4), mar = c(3.1, 3.1, 2.4, 0.8), mgp = c(1.8, 0.6, 0))
    for (j in seq_len(ncol(z))) {
      qqnorm(z[, j], main = paste0(colnames(z)[j], ", n=", nm),
             pch = 16, cex = 0.45, col = "#2C7FB8")
      abline(0, 1, col = "#D7301F", lwd = 1.5)
      qqline(z[, j], col = "grey45", lwd = 1, lty = 2)
    }
    for (j in seq_len(16L - ncol(z))) plot.new()
  }
}

plot_joint_qq <- function(whitened_by_n, file) {
  grDevices::pdf(file, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- ncol(whitened_by_n[[1]])
  old <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  on.exit(par(old), add = TRUE)
  for (nm in names(whitened_by_n)) {
    d2 <- sort(rowSums(whitened_by_n[[nm]]^2))
    theo <- qchisq(ppoints(length(d2)), df = p)
    plot(theo, d2, pch = 16, cex = 0.55, col = "#238B45",
         xlab = expression(chi[p]^2 ~ quantiles),
         ylab = "squared whitened errors", main = paste0("n = ", nm))
    abline(0, 1, col = "#D7301F", lwd = 1.5)
  }
}

run_scenario1c_normality <- function(
    B = 200L, n_values = c(1000L, 2000L, 5000L), seed = NULL,
    setting = "s1c",
    ncores = 1L, nstart = 3L,
    bandwidth_constant = 0.5,
    beta = NULL,
    lambda_constant = 0.10, lambda_power = 0.25, gamma = 3,
    max_iter = 150L, maxeval = 600L,
    information_grid_n = 200001L,
    kde_grid_n_max = 131072L,
    kde_points_per_bandwidth = 10,
    include_population_target = TRUE,
    population_target_maxeval = 2000L,
    population_kernel_nodes = 401L,
    newton_hessian_step = 1e-4,
    local_trust_radius = 8,
    score_ridge_relative = 1e-3,
    score_max_rootn_step = 6,
    score_bootstrap_B = 0L,
    score_bootstrap_min_valid_fraction = 0.7,
    out_dir = file.path(script_dir, "normality_scenario1c_diagnostic_v2")) {
  setting <- tolower(setting)
  if (!setting %in% names(scenario1_settings)) {
    stop("setting must be one of: ", paste(names(scenario1_settings),
                                            collapse = ", "))
  }
  truth <- scenario1_settings[[setting]]
  if (is.null(seed)) seed <- if (setting == "s1b") 211L else 307L
  a_bar <- min(pmin(truth$alpha, 1))
  beta_lower <- 1 / (4 * a_bar)
  if (is.null(beta)) {
    beta <- if (1 / 3 > beta_lower) 1 / 3 else (beta_lower + 0.5) / 2
  }
  if (!(beta > beta_lower && beta < 0.5)) {
    stop("beta must lie in (", signif(beta_lower, 4), ", 0.5) for Scenario 1(c).")
  }
  if (!is.finite(newton_hessian_step) || newton_hessian_step <= 0) {
    stop("newton_hessian_step must be a positive finite number")
  }
  if (!is.finite(local_trust_radius) || local_trust_radius <= 0) {
    stop("local_trust_radius must be a positive finite number")
  }
  if (!is.finite(score_ridge_relative) || score_ridge_relative <= 0 ||
      score_ridge_relative >= 1) {
    stop("score_ridge_relative must lie in (0,1)")
  }
  if (!is.finite(score_max_rootn_step) || score_max_rootn_step <= 0) {
    stop("score_max_rootn_step must be positive")
  }
  score_bootstrap_B <- as.integer(score_bootstrap_B)
  if (score_bootstrap_B < 0L) stop("score_bootstrap_B cannot be negative")
  if (!is.finite(score_bootstrap_min_valid_fraction) ||
      score_bootstrap_min_valid_fraction <= 0 ||
      score_bootstrap_min_valid_fraction > 1) {
    stop("score_bootstrap_min_valid_fraction must lie in (0,1]")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat("Scenario 1 active-submodel normality simulation; setting = ",
      setting, "\n", sep = "")
  cat("B =", B, "; n =", paste(n_values, collapse = ", "), "\n")
  cat("h_n = ", bandwidth_constant, " * n^(-", beta, ")",
      "; admissible beta interval = (", beta_lower, ", 0.5)\n", sep = "")
  cat("Oracle local-branch trust radius =", local_trust_radius,
      "in root-n Fisher-whitened coordinates\n")
  cat("KDE-score one-step: natural chart; ridge =", score_ridge_relative,
      "; max root-n Fisher step =", score_max_rootn_step,
      "; bootstrap draws =", score_bootstrap_B, "\n")
  cat("Computing the theoretical information matrix ...\n")
  J <- information_matrix(truth, grid_n = information_grid_n)
  Sigma <- safe_solve(J)
  chol_J <- chol(J)
  theta0 <- pack_theta(truth)
  captured_mass <- attr(J, "captured_mass")
  information_eigenvalues <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
  information_diagnostics <- data.frame(
    captured_mass = captured_mass,
    minimum_eigenvalue = min(information_eigenvalues),
    maximum_eigenvalue = max(information_eigenvalues),
    spectral_condition_number = max(information_eigenvalues) /
      min(information_eigenvalues)
  )
  cat("Information-grid captured mass:", format(captured_mass, digits = 10), "\n")

  population_target_table <- data.frame()
  if (isTRUE(include_population_target)) {
    cat("Computing deterministic population-smoothed targets ...\n")
    population_target_table <- do.call(rbind, lapply(
      as.integer(n_values), function(n_target) {
        cat("  population target n =", n_target, "\n")
        compute_population_smoothed_target(
          n = n_target, truth = truth, beta = beta,
          bandwidth_constant = bandwidth_constant, Sigma = Sigma,
          maxeval = population_target_maxeval,
          kernel_nodes = population_kernel_nodes,
          grid_n_max = kde_grid_n_max,
          points_per_bandwidth = kde_points_per_bandwidth,
          local_trust_radius = local_trust_radius
        )
      }
    ))
    write.csv(
      population_target_table,
      file.path(out_dir, "population_smoothed_target.csv"),
      row.names = FALSE
    )
  }

  raw_by_n <- list()
  standardized_by_n <- list()
  whitened_by_n <- list()
  linear_raw_by_n <- list()
  linear_standardized_by_n <- list()
  linear_whitened_by_n <- list()
  kde_linear_raw_by_n <- list()
  kde_linear_standardized_by_n <- list()
  kde_linear_whitened_by_n <- list()
  fisher_newton_raw_by_n <- list()
  fisher_newton_standardized_by_n <- list()
  observed_newton_raw_by_n <- list()
  observed_newton_standardized_by_n <- list()
  score_one_step_raw_by_n <- list()
  score_one_step_standardized_by_n <- list()
  score_one_step_whitened_by_n <- list()
  score_bootstrap_full_raw_by_n <- list()
  score_bootstrap_full_standardized_by_n <- list()
  score_bootstrap_full_whitened_by_n <- list()
  score_bootstrap_shrunk_raw_by_n <- list()
  score_bootstrap_shrunk_standardized_by_n <- list()
  score_bootstrap_shrunk_whitened_by_n <- list()
  score_one_step_diagnostics <- list()
  score_bootstrap_diagnostics <- list()
  remainder_by_n <- list()
  newton_diagnostics <- list()
  fit_diagnostics <- list()
  failure_diagnostics <- list()
  all_results <- list()

  cl <- NULL
  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    core_file_worker <- core_file
    parallel::clusterExport(cl, "core_file_worker", envir = environment())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(source(core_file_worker))
      NULL
    })
  }

  for (n in as.integer(n_values)) {
    cat("\nRunning n =", n, "...\n")
    args <- list(n = n, truth = truth, beta = beta,
                 bandwidth_constant = bandwidth_constant,
                 lambda_constant = lambda_constant,
                 lambda_power = lambda_power, gamma = gamma,
                 nstart = nstart, max_iter = max_iter, maxeval = maxeval,
                 Sigma = Sigma, J = J,
                 kde_grid_n_max = kde_grid_n_max,
                 kde_points_per_bandwidth = kde_points_per_bandwidth,
                 newton_hessian_step = newton_hessian_step,
                 local_trust_radius = local_trust_radius,
                 score_ridge_relative = score_ridge_relative,
                 score_max_rootn_step = score_max_rootn_step,
                 score_bootstrap_B = score_bootstrap_B,
                 score_bootstrap_min_valid_fraction =
                   score_bootstrap_min_valid_fraction)
    if (is.null(cl)) {
      set.seed(seed + n)
      res <- lapply(seq_len(B), function(b) {
        if (b %% max(1L, B %/% 10L) == 0L) cat("  replication", b, "of", B, "\n")
        do.call(simulate_one_normality, c(list(b = b), args))
      })
    } else {
      parallel::clusterSetRNGStream(cl, seed + n)
      parallel::clusterExport(
        cl,
        c("simulate_one_normality", "fit_active_mhd", "make_theory_kde",
          "make_local_starts", "align_by_location", "pack_theta",
          "parameter_names", "softmax", "next_power_of_two",
          "working_from_natural", "natural_from_working",
          "hellinger_objective_gradient", "direct_refine_active_mhd",
          "direct_refine_active_mhd_wide",
          "natural_jacobian_from_working",
          "finite_difference_gradient_hessian",
          "compute_kde_newton_diagnostics", "safe_solve",
          "aepd_component_scores", "mixture_score_matrix",
          "kde_score_one_step", "unpack_theta_natural",
          "kde_score_bootstrap_centering", "args"),
        envir = environment()
      )
      res <- parallel::parLapply(cl, seq_len(B), function(b) {
        do.call(simulate_one_normality, c(list(b = b), args))
      })
    }
    ok <- vapply(res, function(z) isTRUE(z$ok), logical(1))
    if (sum(ok) < max(2L, ceiling(0.2 * B))) {
      stop("Only ", sum(ok), "/", B, " fits succeeded at n=", n,
           ". Inspect initialization or optimizer settings.")
    }
    failure_reason <- vapply(res[!ok], function(z) {
      if (is.null(z$message)) "unknown" else z$message
    }, character(1))
    if (length(failure_reason)) {
      failure_diagnostics[[as.character(n)]] <- data.frame(
        n = n, reason = names(table(failure_reason)),
        count = as.integer(table(failure_reason)), row.names = NULL
      )
    }
    good <- res[ok]
    linear_rootn_error_all <- do.call(
      rbind, lapply(res, `[[`, "linear_rootn_error")
    )
    kde_linear_rootn_error_all <- do.call(
      rbind, lapply(res, `[[`, "kde_linear_rootn_error")
    )
    fisher_newton_all <- do.call(
      rbind, lapply(res, `[[`, "fisher_rootn_error")
    )
    observed_newton_all <- do.call(
      rbind, lapply(res, `[[`, "observed_rootn_error")
    )
    colnames(linear_rootn_error_all) <-
      colnames(kde_linear_rootn_error_all) <- names(theta0)
    colnames(fisher_newton_all) <- colnames(observed_newton_all) <-
      names(theta0)
    linear_marginal_z <- sweep(
      linear_rootn_error_all, 2, sqrt(diag(Sigma)), "/"
    )
    kde_linear_marginal_z <- sweep(
      kde_linear_rootn_error_all, 2, sqrt(diag(Sigma)), "/"
    )
    fisher_newton_z <- sweep(
      fisher_newton_all, 2, sqrt(diag(Sigma)), "/"
    )
    observed_newton_z <- sweep(
      observed_newton_all, 2, sqrt(diag(Sigma)), "/"
    )
    linear_whitened <- linear_rootn_error_all %*% t(chol_J)
    kde_linear_whitened <- kde_linear_rootn_error_all %*% t(chol_J)
    colnames(linear_marginal_z) <- colnames(linear_whitened) <-
      colnames(kde_linear_marginal_z) <- colnames(kde_linear_whitened) <-
      names(theta0)
    estimates <- do.call(rbind, lapply(good, `[[`, "theta"))
    root_n_error <- sqrt(n) * sweep(estimates, 2, theta0, "-")
    marginal_z <- sweep(root_n_error, 2, sqrt(diag(Sigma)), "/")
    whitened <- root_n_error %*% t(chol_J)
    colnames(root_n_error) <- colnames(marginal_z) <- colnames(whitened) <- names(theta0)
    score_one_step_estimates <- do.call(
      rbind, lapply(good, `[[`, "kde_score_one_step_theta")
    )
    colnames(score_one_step_estimates) <- names(theta0)
    score_one_step_rootn_error <-
      sqrt(n) * sweep(score_one_step_estimates, 2, theta0, "-")
    score_one_step_z <- sweep(
      score_one_step_rootn_error, 2, sqrt(diag(Sigma)), "/"
    )
    score_one_step_whitened <- score_one_step_rootn_error %*% t(chol_J)
    colnames(score_one_step_rootn_error) <-
      colnames(score_one_step_z) <- colnames(score_one_step_whitened) <-
      names(theta0)
    if (score_bootstrap_B > 0L) {
      score_bootstrap_full_estimates <- do.call(
        rbind, lapply(good, `[[`, "kde_score_bootstrap_full_theta")
      )
      score_bootstrap_shrunk_estimates <- do.call(
        rbind, lapply(good, `[[`, "kde_score_bootstrap_shrunk_theta")
      )
      colnames(score_bootstrap_full_estimates) <-
        colnames(score_bootstrap_shrunk_estimates) <- names(theta0)
      score_bootstrap_full_rootn_error <- sqrt(n) * sweep(
        score_bootstrap_full_estimates, 2, theta0, "-"
      )
      score_bootstrap_shrunk_rootn_error <- sqrt(n) * sweep(
        score_bootstrap_shrunk_estimates, 2, theta0, "-"
      )
      score_bootstrap_full_z <- sweep(
        score_bootstrap_full_rootn_error, 2, sqrt(diag(Sigma)), "/"
      )
      score_bootstrap_shrunk_z <- sweep(
        score_bootstrap_shrunk_rootn_error, 2, sqrt(diag(Sigma)), "/"
      )
      score_bootstrap_full_whitened <-
        score_bootstrap_full_rootn_error %*% t(chol_J)
      score_bootstrap_shrunk_whitened <-
        score_bootstrap_shrunk_rootn_error %*% t(chol_J)
      colnames(score_bootstrap_full_rootn_error) <-
        colnames(score_bootstrap_shrunk_rootn_error) <-
        colnames(score_bootstrap_full_z) <-
        colnames(score_bootstrap_shrunk_z) <-
        colnames(score_bootstrap_full_whitened) <-
        colnames(score_bootstrap_shrunk_whitened) <- names(theta0)
    }

    key <- as.character(n)
    raw_by_n[[key]] <- root_n_error
    standardized_by_n[[key]] <- marginal_z
    whitened_by_n[[key]] <- whitened
    linear_raw_by_n[[key]] <- linear_rootn_error_all
    linear_standardized_by_n[[key]] <- linear_marginal_z
    linear_whitened_by_n[[key]] <- linear_whitened
    kde_linear_raw_by_n[[key]] <- kde_linear_rootn_error_all
    kde_linear_standardized_by_n[[key]] <- kde_linear_marginal_z
    kde_linear_whitened_by_n[[key]] <- kde_linear_whitened
    fisher_newton_raw_by_n[[key]] <- fisher_newton_all
    fisher_newton_standardized_by_n[[key]] <- fisher_newton_z
    observed_newton_raw_by_n[[key]] <- observed_newton_all
    observed_newton_standardized_by_n[[key]] <- observed_newton_z
    score_one_step_raw_by_n[[key]] <- score_one_step_rootn_error
    score_one_step_standardized_by_n[[key]] <- score_one_step_z
    score_one_step_whitened_by_n[[key]] <- score_one_step_whitened
    if (score_bootstrap_B > 0L) {
      score_bootstrap_full_raw_by_n[[key]] <-
        score_bootstrap_full_rootn_error
      score_bootstrap_full_standardized_by_n[[key]] <-
        score_bootstrap_full_z
      score_bootstrap_full_whitened_by_n[[key]] <-
        score_bootstrap_full_whitened
      score_bootstrap_shrunk_raw_by_n[[key]] <-
        score_bootstrap_shrunk_rootn_error
      score_bootstrap_shrunk_standardized_by_n[[key]] <-
        score_bootstrap_shrunk_z
      score_bootstrap_shrunk_whitened_by_n[[key]] <-
        score_bootstrap_shrunk_whitened
    }
    score_one_step_diagnostics[[key]] <- data.frame(
      n = n, outer_index = which(ok),
      moment_max_abs = vapply(
        good, `[[`, numeric(1), "kde_score_moment_max_abs"
      ),
      moment_l2 = vapply(good, `[[`, numeric(1), "kde_score_moment_l2"),
      rootn_fisher_step = vapply(
        good, `[[`, numeric(1), "kde_score_rootn_fisher_step"
      ),
      information_min_eigenvalue = vapply(
        good, `[[`, numeric(1), "kde_score_information_min_eigenvalue"
      ),
      information_condition_number = vapply(
        good, `[[`, numeric(1), "kde_score_information_condition_number"
      ),
      information_ridge_used = vapply(
        good, `[[`, logical(1), "kde_score_information_ridge_used"
      ),
      model_mass_on_grid = vapply(
        good, `[[`, numeric(1), "kde_score_model_mass_on_grid"
      ),
      step_clipped = vapply(
        good, `[[`, logical(1), "kde_score_step_clipped"
      ),
      step_scale = vapply(good, `[[`, numeric(1), "kde_score_step_scale"),
      step_halvings = vapply(
        good, `[[`, integer(1), "kde_score_step_halvings"
      ),
      bootstrap_valid = vapply(
        good, `[[`, integer(1), "kde_score_bootstrap_valid"
      ),
      bootstrap_full_scale = vapply(
        good, `[[`, numeric(1), "kde_score_bootstrap_full_scale"
      ),
      bootstrap_shrunk_scale = vapply(
        good, `[[`, numeric(1), "kde_score_bootstrap_shrunk_scale"
      )
    )
    if (score_bootstrap_B > 0L) {
      bootstrap_bias <- do.call(
        rbind, lapply(good, `[[`, "kde_score_bootstrap_bias_rootn")
      )
      bootstrap_bias_mcse <- do.call(
        rbind, lapply(good, `[[`, "kde_score_bootstrap_bias_mcse_rootn")
      )
      bootstrap_shrinkage <- do.call(
        rbind, lapply(good, `[[`, "kde_score_bootstrap_shrinkage")
      )
      score_bootstrap_diagnostics[[key]] <- data.frame(
        n = n,
        outer_index = rep(which(ok), each = length(theta0)),
        parameter = rep(names(theta0), times = length(good)),
        bias_rootn = as.vector(t(bootstrap_bias)),
        bias_mcse_rootn = as.vector(t(bootstrap_bias_mcse)),
        shrinkage = as.vector(t(bootstrap_shrinkage)),
        valid_bootstrap_draws = rep(vapply(
          good, `[[`, integer(1), "kde_score_bootstrap_valid"
        ), each = length(theta0)),
        full_correction_scale = rep(vapply(
          good, `[[`, numeric(1), "kde_score_bootstrap_full_scale"
        ), each = length(theta0)),
        shrunk_correction_scale = rep(vapply(
          good, `[[`, numeric(1), "kde_score_bootstrap_shrunk_scale"
        ), each = length(theta0))
      )
    }
    good_index <- which(ok)
    remainder_by_n[[key]] <- root_n_error -
      kde_linear_rootn_error_all[good_index, , drop = FALSE]
    newton_diagnostics[[key]] <- data.frame(
      n = n,
      objective_at_truth = vapply(res, `[[`, numeric(1),
                                  "objective_at_truth"),
      gradient_max_abs_at_truth = vapply(
        res, `[[`, numeric(1), "gradient_max_abs_at_truth"
      ),
      gradient_l2_at_truth = vapply(
        res, `[[`, numeric(1), "gradient_l2_at_truth"
      ),
      observed_hessian_min_eigenvalue = vapply(
        res, `[[`, numeric(1), "observed_hessian_min_eigenvalue"
      ),
      observed_hessian_max_eigenvalue = vapply(
        res, `[[`, numeric(1), "observed_hessian_max_eigenvalue"
      ),
      observed_hessian_condition_number = vapply(
        res, `[[`, numeric(1), "observed_hessian_condition_number"
      ),
      observed_hessian_positive_definite = vapply(
        res, `[[`, logical(1), "observed_hessian_positive_definite"
      ),
      hessian_location_step_over_grid = vapply(
        res, `[[`, numeric(1), "hessian_location_step_over_grid"
      )
    )
    fit_diagnostics[[key]] <- data.frame(
      n = n, requested = B, successful = sum(ok), failed = sum(!ok),
      flat_tail_rate = mean(vapply(good, `[[`, logical(1), "flat_tail")),
      direct_success_rate = mean(vapply(good, function(z) {
        is.finite(z$optimizer_status) && z$optimizer_status > 0
      }, logical(1))),
      median_gradient_max_abs = median(vapply(
        good, `[[`, numeric(1), "gradient_max_abs"
      )),
      median_local_gradient_max_abs = median(vapply(
        good, `[[`, numeric(1), "local_gradient_max_abs"
      )),
      median_trust_coordinate_max_abs = median(vapply(
        good, `[[`, numeric(1), "trust_coordinate_max_abs"
      )),
      max_trust_coordinate_max_abs = max(vapply(
        good, `[[`, numeric(1), "trust_coordinate_max_abs"
      )),
      median_distance_from_fisher_l2 = median(vapply(
        good, `[[`, numeric(1), "distance_from_fisher_l2"
      )),
      mean_local_candidates_valid = mean(vapply(
        good, `[[`, numeric(1), "local_candidates_valid"
      )),
      median_objective_decrease_from_truth = median(vapply(
        good, `[[`, numeric(1), "objective_decrease_from_truth"
      )),
      local_hessian_pd_rate = mean(vapply(
        good, `[[`, logical(1), "hessian_positive_definite"
      )),
      median_local_hessian_min_eigenvalue = median(vapply(
        good, `[[`, numeric(1), "hessian_min_eigenvalue"
      )),
      alpha_boundary_rate = mean(vapply(
        good, `[[`, logical(1), "boundary_alpha"
      )),
      median_grid_n = median(vapply(good, `[[`, numeric(1), "grid_n")),
      max_grid_n = max(vapply(good, `[[`, numeric(1), "grid_n")),
      mean_grid_step_over_h = mean(vapply(
        good, `[[`, numeric(1), "grid_step_over_h"
      )),
      max_grid_step_over_h = max(vapply(
        good, `[[`, numeric(1), "grid_step_over_h"
      )),
      grid_cap_hit_rate = mean(vapply(
        good, `[[`, logical(1), "grid_cap_hit"
      )),
      mean_kde_mass_before_normalizing = mean(vapply(
        good, `[[`, numeric(1), "kde_mass_before_normalizing"
      )),
      mean_H2 = mean(vapply(good, `[[`, numeric(1), "H2")),
      newton_hessian_pd_rate = mean(vapply(
        res, `[[`, logical(1), "observed_hessian_positive_definite"
      )),
      median_newton_hessian_condition_number = median(vapply(
        res, `[[`, numeric(1), "observed_hessian_condition_number"
      ), na.rm = TRUE),
      bandwidth = good[[1]]$h, lambda = good[[1]]$lambda,
      local_trust_radius = local_trust_radius
    )
    all_results[[key]] <- res
  }

  summary_table <- do.call(rbind, Map(normal_summary, standardized_by_n,
                                      as.integer(names(standardized_by_n)),
                                      MoreArgs = list(names_z = names(theta0))))
  diagnostics_table <- do.call(rbind, fit_diagnostics)
  success_rate <- setNames(diagnostics_table$successful /
                             diagnostics_table$requested,
                           diagnostics_table$n)
  summary_table$success_rate <- success_rate[as.character(summary_table$n)]
  summary_table$coverage_95_unconditional <-
    summary_table$coverage_95 * summary_table$success_rate
  linear_summary_table <- do.call(rbind, Map(
    normal_summary, linear_standardized_by_n,
    as.integer(names(linear_standardized_by_n)),
    MoreArgs = list(names_z = names(theta0))
  ))
  linear_summary_table$success_rate <- 1
  linear_summary_table$coverage_95_unconditional <-
    linear_summary_table$coverage_95
  kde_linear_summary_table <- do.call(rbind, lapply(
    names(kde_linear_standardized_by_n), function(nm) {
      z <- kde_linear_standardized_by_n[[nm]]
      z <- z[complete.cases(z), , drop = FALSE]
      normal_summary(z, as.integer(nm), names(theta0))
    }
  ))
  kde_linear_summary_table$success_rate <- 1
  kde_linear_summary_table$coverage_95_unconditional <-
    kde_linear_summary_table$coverage_95
  fisher_newton_summary_table <- do.call(rbind, lapply(
    names(fisher_newton_standardized_by_n), function(nm) {
      z <- fisher_newton_standardized_by_n[[nm]]
      z <- z[complete.cases(z), , drop = FALSE]
      normal_summary(z, as.integer(nm), names(theta0))
    }
  ))
  observed_newton_summary_table <- do.call(rbind, lapply(
    names(observed_newton_standardized_by_n), function(nm) {
      z <- observed_newton_standardized_by_n[[nm]]
      z <- z[complete.cases(z), , drop = FALSE]
      normal_summary(z, as.integer(nm), names(theta0))
    }
  ))
  score_one_step_summary_table <- do.call(rbind, lapply(
    names(score_one_step_standardized_by_n), function(nm) {
      z <- score_one_step_standardized_by_n[[nm]]
      z <- z[complete.cases(z), , drop = FALSE]
      normal_summary(z, as.integer(nm), names(theta0))
    }
  ))
  score_one_step_summary_table$success_rate <-
    success_rate[as.character(score_one_step_summary_table$n)]
  score_one_step_summary_table$coverage_95_unconditional <-
    score_one_step_summary_table$coverage_95 *
    score_one_step_summary_table$success_rate
  bootstrap_summary <- function(zlist) {
    tab <- do.call(rbind, lapply(names(zlist), function(nm) {
      z <- zlist[[nm]]
      z <- z[complete.cases(z), , drop = FALSE]
      normal_summary(z, as.integer(nm), names(theta0))
    }))
    tab$success_rate <- success_rate[as.character(tab$n)]
    tab$coverage_95_unconditional <- tab$coverage_95 * tab$success_rate
    tab
  }
  score_bootstrap_full_summary_table <-
    score_bootstrap_shrunk_summary_table <- data.frame()
  if (score_bootstrap_B > 0L) {
    score_bootstrap_full_summary_table <- bootstrap_summary(
      score_bootstrap_full_standardized_by_n
    )
    score_bootstrap_shrunk_summary_table <- bootstrap_summary(
      score_bootstrap_shrunk_standardized_by_n
    )
  }
  summary_table <- add_standard_deviation_comparison(
    summary_table, Sigma, names(theta0)
  )
  linear_summary_table <- add_standard_deviation_comparison(
    linear_summary_table, Sigma, names(theta0)
  )
  kde_linear_summary_table <- add_standard_deviation_comparison(
    kde_linear_summary_table, Sigma, names(theta0)
  )
  fisher_newton_summary_table <- add_standard_deviation_comparison(
    fisher_newton_summary_table, Sigma, names(theta0)
  )
  observed_newton_summary_table <- add_standard_deviation_comparison(
    observed_newton_summary_table, Sigma, names(theta0)
  )
  score_one_step_summary_table <- add_standard_deviation_comparison(
    score_one_step_summary_table, Sigma, names(theta0)
  )
  if (score_bootstrap_B > 0L) {
    score_bootstrap_full_summary_table <- add_standard_deviation_comparison(
      score_bootstrap_full_summary_table, Sigma, names(theta0)
    )
    score_bootstrap_shrunk_summary_table <- add_standard_deviation_comparison(
      score_bootstrap_shrunk_summary_table, Sigma, names(theta0)
    )
  }
  standard_deviation_comparison_table <- summary_table[, c(
    "n", "parameter", "replications", "success_rate",
    "empirical_sd", "theoretical_sd",
    "empirical_to_theoretical_sd_ratio",
    "coverage_95", "coverage_95_mcse",
    "empirical_root_n_sd", "theoretical_root_n_sd", "standardized_sd"
  )]
  marginal_sd_coverage_table <- summary_table[, c(
    "n", "parameter", "replications", "success_rate",
    "empirical_sd", "theoretical_sd", "coverage_95"
  )]
  remainder_table <- do.call(rbind, lapply(names(remainder_by_n), function(nm) {
    rem <- remainder_by_n[[nm]]
    actual <- raw_by_n[[nm]]
    idx <- which(vapply(all_results[[nm]], function(z) isTRUE(z$ok), logical(1)))
    linear <- kde_linear_raw_by_n[[nm]][idx, , drop = FALSE]
    data.frame(
      n = as.integer(nm), parameter = colnames(rem),
      mean_remainder = colMeans(rem),
      sd_remainder = apply(rem, 2, sd),
      rmse_remainder = sqrt(colMeans(rem^2)),
      relative_rmse = sqrt(colMeans(rem^2)) /
        pmax(sqrt(colMeans(linear^2)), 1e-12),
      estimator_linear_correlation = vapply(seq_len(ncol(rem)), function(j) {
        cor(actual[, j], linear[, j])
      }, numeric(1))
    )
  }))
  score_one_step_remainder_table <- do.call(rbind, lapply(
    names(score_one_step_raw_by_n), function(nm) {
      actual <- score_one_step_raw_by_n[[nm]]
      idx <- which(vapply(
        all_results[[nm]], function(z) isTRUE(z$ok), logical(1)
      ))
      linear <- kde_linear_raw_by_n[[nm]][idx, , drop = FALSE]
      rem <- actual - linear
      data.frame(
        n = as.integer(nm), parameter = colnames(rem),
        mean_remainder = colMeans(rem, na.rm = TRUE),
        sd_remainder = apply(rem, 2, sd, na.rm = TRUE),
        rmse_remainder = sqrt(colMeans(rem^2, na.rm = TRUE)),
        relative_rmse = sqrt(colMeans(rem^2, na.rm = TRUE)) /
          pmax(sqrt(colMeans(linear^2, na.rm = TRUE)), 1e-12),
        one_step_linear_correlation = vapply(
          seq_len(ncol(rem)), function(j) {
            cor(actual[, j], linear[, j], use = "complete.obs")
          }, numeric(1)
        )
      )
    }
  ))
  newton_comparison_table <- do.call(rbind, lapply(names(raw_by_n), function(nm) {
    se <- sqrt(diag(Sigma))
    linear_z <- sweep(linear_raw_by_n[[nm]], 2, se, "/")
    kde_linear_z <- sweep(kde_linear_raw_by_n[[nm]], 2, se, "/")
    fisher_z <- sweep(fisher_newton_raw_by_n[[nm]], 2, se, "/")
    observed_z <- sweep(observed_newton_raw_by_n[[nm]], 2, se, "/")
    good_index <- which(vapply(
      all_results[[nm]], function(z) isTRUE(z$ok), logical(1)
    ))
    estimator_z <- standardized_by_n[[nm]]
    fisher_good_z <- fisher_z[good_index, , drop = FALSE]
    observed_good_z <- observed_z[good_index, , drop = FALSE]
    data.frame(
      n = as.integer(nm), parameter = names(theta0),
      fisher_mean_z = colMeans(fisher_z, na.rm = TRUE),
      fisher_sd_z = apply(fisher_z, 2, sd, na.rm = TRUE),
      kde_linear_mean_z = colMeans(kde_linear_z, na.rm = TRUE),
      kde_linear_sd_z = apply(kde_linear_z, 2, sd, na.rm = TRUE),
      kde_raw_score_rmse_z = sqrt(colMeans(
        (kde_linear_z - linear_z)^2, na.rm = TRUE
      )),
      kde_raw_score_correlation = vapply(seq_along(theta0), function(j) {
        cor(kde_linear_z[, j], linear_z[, j], use = "complete.obs")
      }, numeric(1)),
      fisher_kde_rmse_z = sqrt(colMeans(
        (fisher_z - kde_linear_z)^2, na.rm = TRUE
      )),
      fisher_kde_correlation = vapply(seq_along(theta0), function(j) {
        cor(fisher_z[, j], kde_linear_z[, j], use = "complete.obs")
      }, numeric(1)),
      fisher_linear_rmse_z = sqrt(colMeans(
        (fisher_z - linear_z)^2, na.rm = TRUE
      )),
      fisher_linear_correlation = vapply(seq_along(theta0), function(j) {
        cor(fisher_z[, j], linear_z[, j], use = "complete.obs")
      }, numeric(1)),
      observed_mean_z = colMeans(observed_z, na.rm = TRUE),
      observed_sd_z = apply(observed_z, 2, sd, na.rm = TRUE),
      observed_fisher_rmse_z = sqrt(colMeans(
        (observed_z - fisher_z)^2, na.rm = TRUE
      )),
      observed_fisher_correlation = vapply(seq_along(theta0), function(j) {
        cor(observed_z[, j], fisher_z[, j], use = "complete.obs")
      }, numeric(1)),
      estimator_fisher_mean_difference_z = colMeans(
        estimator_z - fisher_good_z, na.rm = TRUE
      ),
      estimator_fisher_rmse_z = sqrt(colMeans(
        (estimator_z - fisher_good_z)^2, na.rm = TRUE
      )),
      estimator_fisher_correlation = vapply(seq_along(theta0), function(j) {
        cor(estimator_z[, j], fisher_good_z[, j], use = "complete.obs")
      }, numeric(1)),
      estimator_observed_rmse_z = sqrt(colMeans(
        (estimator_z - observed_good_z)^2, na.rm = TRUE
      )),
      estimator_observed_correlation = vapply(seq_along(theta0), function(j) {
        cor(estimator_z[, j], observed_good_z[, j], use = "complete.obs")
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  newton_diagnostics_table <- do.call(rbind, newton_diagnostics)
  failure_table <- if (length(failure_diagnostics)) {
    do.call(rbind, failure_diagnostics)
  } else {
    data.frame(n = integer(), reason = character(), count = integer())
  }
  covariance_table <- do.call(rbind, lapply(names(raw_by_n), function(nm) {
    empirical <- stats::cov(raw_by_n[[nm]])
    squared_whitened_error <- rowSums(whitened_by_n[[nm]]^2)
    joint_covered_95 <- squared_whitened_error <=
      qchisq(0.95, df = ncol(whitened_by_n[[nm]]))
    data.frame(
      n = as.integer(nm),
      relative_frobenius_error = sqrt(sum((empirical - Sigma)^2)) /
        sqrt(sum(Sigma^2)),
      joint_qq_correlation = cor(
        sort(squared_whitened_error),
        qchisq(ppoints(nrow(whitened_by_n[[nm]])), df = ncol(whitened_by_n[[nm]]))
      ),
      joint_coverage_95 = mean(joint_covered_95),
      success_rate = success_rate[nm],
      joint_coverage_95_unconditional =
        mean(joint_covered_95) * success_rate[nm],
      joint_coverage_95_mcse = sqrt(
        mean(joint_covered_95) * (1 - mean(joint_covered_95)) /
          length(joint_covered_95)
      )
    )
  }))
  linear_joint_table <- do.call(rbind, lapply(names(linear_raw_by_n), function(nm) {
    empirical <- stats::cov(linear_raw_by_n[[nm]])
    squared_whitened_error <- rowSums(linear_whitened_by_n[[nm]]^2)
    hit <- squared_whitened_error <= qchisq(
      0.95, df = ncol(linear_whitened_by_n[[nm]])
    )
    data.frame(
      n = as.integer(nm),
      relative_frobenius_error = sqrt(sum((empirical - Sigma)^2)) /
        sqrt(sum(Sigma^2)),
      joint_qq_correlation = cor(
        sort(squared_whitened_error),
        qchisq(ppoints(length(squared_whitened_error)),
               df = ncol(linear_whitened_by_n[[nm]]))
      ),
      joint_coverage_95 = mean(hit),
      joint_coverage_95_mcse = sqrt(mean(hit) * (1 - mean(hit)) / length(hit))
    )
  }))
  kde_linear_joint_table <- do.call(rbind, lapply(
    names(kde_linear_raw_by_n), function(nm) {
      raw <- kde_linear_raw_by_n[[nm]]
      white <- kde_linear_whitened_by_n[[nm]]
      keep <- complete.cases(raw, white)
      raw <- raw[keep, , drop = FALSE]
      white <- white[keep, , drop = FALSE]
      empirical <- stats::cov(raw)
      squared_whitened_error <- rowSums(white^2)
      hit <- squared_whitened_error <= qchisq(
        0.95, df = ncol(white)
      )
      data.frame(
        n = as.integer(nm),
        replications = nrow(raw),
        relative_frobenius_error = sqrt(sum((empirical - Sigma)^2)) /
          sqrt(sum(Sigma^2)),
        joint_qq_correlation = cor(
          sort(squared_whitened_error),
          qchisq(ppoints(length(squared_whitened_error)),
                 df = ncol(white))
        ),
        joint_coverage_95 = mean(hit),
        joint_coverage_95_mcse =
          sqrt(mean(hit) * (1 - mean(hit)) / length(hit))
      )
    }
  ))
  score_one_step_joint_table <- do.call(rbind, lapply(
    names(score_one_step_raw_by_n), function(nm) {
      raw <- score_one_step_raw_by_n[[nm]]
      white <- score_one_step_whitened_by_n[[nm]]
      keep <- complete.cases(raw, white)
      raw <- raw[keep, , drop = FALSE]
      white <- white[keep, , drop = FALSE]
      empirical <- stats::cov(raw)
      squared_whitened_error <- rowSums(white^2)
      hit <- squared_whitened_error <= qchisq(0.95, df = ncol(white))
      data.frame(
        n = as.integer(nm), replications = nrow(raw),
        relative_frobenius_error = sqrt(sum((empirical - Sigma)^2)) /
          sqrt(sum(Sigma^2)),
        joint_qq_correlation = cor(
          sort(squared_whitened_error),
          qchisq(ppoints(length(squared_whitened_error)), df = ncol(white))
        ),
        joint_coverage_95 = mean(hit),
        success_rate = success_rate[nm],
        joint_coverage_95_unconditional = mean(hit) * success_rate[nm],
        joint_coverage_95_mcse =
          sqrt(mean(hit) * (1 - mean(hit)) / length(hit))
      )
    }
  ))
  bootstrap_joint_summary <- function(raw_by_n, whitened_by_n) {
    do.call(rbind, lapply(names(raw_by_n), function(nm) {
      raw <- raw_by_n[[nm]]
      white <- whitened_by_n[[nm]]
      keep <- complete.cases(raw, white)
      raw <- raw[keep, , drop = FALSE]
      white <- white[keep, , drop = FALSE]
      if (nrow(raw) < 2L) {
        return(data.frame(
          n = as.integer(nm), replications = nrow(raw),
          relative_frobenius_error = NA_real_,
          joint_qq_correlation = NA_real_, joint_coverage_95 = NA_real_,
          success_rate = success_rate[nm],
          joint_coverage_95_unconditional = NA_real_,
          joint_coverage_95_mcse = NA_real_
        ))
      }
      empirical <- stats::cov(raw)
      squared_whitened_error <- rowSums(white^2)
      hit <- squared_whitened_error <= qchisq(0.95, df = ncol(white))
      data.frame(
        n = as.integer(nm), replications = nrow(raw),
        relative_frobenius_error = sqrt(sum((empirical - Sigma)^2)) /
          sqrt(sum(Sigma^2)),
        joint_qq_correlation = cor(
          sort(squared_whitened_error),
          qchisq(ppoints(length(squared_whitened_error)), df = ncol(white))
        ),
        joint_coverage_95 = mean(hit),
        success_rate = success_rate[nm],
        joint_coverage_95_unconditional = mean(hit) * success_rate[nm],
        joint_coverage_95_mcse =
          sqrt(mean(hit) * (1 - mean(hit)) / length(hit))
      )
    }))
  }
  score_bootstrap_full_joint_table <-
    score_bootstrap_shrunk_joint_table <- data.frame()
  if (score_bootstrap_B > 0L) {
    score_bootstrap_full_joint_table <- bootstrap_joint_summary(
      score_bootstrap_full_raw_by_n, score_bootstrap_full_whitened_by_n
    )
    score_bootstrap_shrunk_joint_table <- bootstrap_joint_summary(
      score_bootstrap_shrunk_raw_by_n, score_bootstrap_shrunk_whitened_by_n
    )
  }
  score_one_step_diagnostics_table <- do.call(
    rbind, score_one_step_diagnostics
  )
  score_bootstrap_diagnostics_table <- data.frame()
  if (score_bootstrap_B > 0L && length(score_bootstrap_diagnostics)) {
    score_bootstrap_diagnostics_table <- do.call(
      rbind, score_bootstrap_diagnostics
    )
  }

  bias_decomposition_table <- data.frame()
  if (nrow(population_target_table)) {
    bias_decomposition_table <- merge(
      summary_table[, c("n", "parameter", "mean", "sd", "coverage_95")],
      population_target_table[, c(
        "n", "parameter", "bandwidth", "theta_true", "theta_h",
        "smoothing_bias", "root_n_smoothing_bias",
        "standardized_smoothing_bias"
      )],
      by = c("n", "parameter"), all.x = TRUE, sort = FALSE
    )
    names(bias_decomposition_table)[
      names(bias_decomposition_table) == "mean"
    ] <- "actual_standardized_mean"
    names(bias_decomposition_table)[
      names(bias_decomposition_table) == "sd"
    ] <- "actual_standardized_sd"
    bias_decomposition_table$residual_standardized_bias <-
      bias_decomposition_table$actual_standardized_mean -
      bias_decomposition_table$standardized_smoothing_bias
  }

  write.csv(summary_table, file.path(out_dir, "marginal_normality_summary.csv"),
            row.names = FALSE)
  write.csv(standard_deviation_comparison_table,
            file.path(out_dir, "standard_deviation_comparison.csv"),
            row.names = FALSE)
  write.csv(marginal_sd_coverage_table,
            file.path(out_dir, "marginal_sd_coverage_summary.csv"),
            row.names = FALSE)
  write.csv(summary_table[, c("n", "parameter", "replications",
                              "success_rate", "coverage_95",
                              "coverage_95_unconditional",
                              "coverage_95_mcse")],
            file.path(out_dir, "marginal_coverage_95.csv"),
            row.names = FALSE)
  write.csv(diagnostics_table, file.path(out_dir, "fit_diagnostics.csv"),
            row.names = FALSE)
  write.csv(failure_table, file.path(out_dir, "failure_diagnostics.csv"),
            row.names = FALSE)
  write.csv(linear_summary_table,
            file.path(out_dir, "linear_benchmark_marginal_summary.csv"),
            row.names = FALSE)
  write.csv(linear_joint_table,
            file.path(out_dir, "linear_benchmark_joint_summary.csv"),
            row.names = FALSE)
  write.csv(kde_linear_summary_table,
            file.path(out_dir, "kde_linear_marginal_summary.csv"),
            row.names = FALSE)
  write.csv(kde_linear_joint_table,
            file.path(out_dir, "kde_linear_joint_summary.csv"),
            row.names = FALSE)
  write.csv(fisher_newton_summary_table,
            file.path(out_dir, "kde_fisher_newton_marginal_summary.csv"),
            row.names = FALSE)
  write.csv(observed_newton_summary_table,
            file.path(out_dir, "kde_observed_newton_marginal_summary.csv"),
            row.names = FALSE)
  write.csv(score_one_step_summary_table,
            file.path(out_dir, "kde_score_one_step_marginal_summary.csv"),
            row.names = FALSE)
  write.csv(score_one_step_joint_table,
            file.path(out_dir, "kde_score_one_step_joint_summary.csv"),
            row.names = FALSE)
  write.csv(score_one_step_diagnostics_table,
            file.path(out_dir, "kde_score_one_step_diagnostics.csv"),
            row.names = FALSE)
  write.csv(score_one_step_remainder_table,
            file.path(out_dir, "kde_score_one_step_remainder.csv"),
            row.names = FALSE)
  if (score_bootstrap_B > 0L) {
    write.csv(score_bootstrap_full_summary_table,
              file.path(out_dir,
                        "kde_score_bootstrap_full_marginal_summary.csv"),
              row.names = FALSE)
    write.csv(score_bootstrap_shrunk_summary_table,
              file.path(out_dir,
                        "kde_score_bootstrap_shrunk_marginal_summary.csv"),
              row.names = FALSE)
    write.csv(score_bootstrap_full_joint_table,
              file.path(out_dir,
                        "kde_score_bootstrap_full_joint_summary.csv"),
              row.names = FALSE)
    write.csv(score_bootstrap_shrunk_joint_table,
              file.path(out_dir,
                        "kde_score_bootstrap_shrunk_joint_summary.csv"),
              row.names = FALSE)
    write.csv(score_bootstrap_diagnostics_table,
              file.path(out_dir, "kde_score_bootstrap_diagnostics.csv"),
              row.names = FALSE)
  }
  write.csv(newton_comparison_table,
            file.path(out_dir, "newton_decomposition_summary.csv"),
            row.names = FALSE)
  write.csv(newton_diagnostics_table,
            file.path(out_dir, "newton_target_diagnostics.csv"),
            row.names = FALSE)
  write.csv(remainder_table,
            file.path(out_dir, "asymptotic_linearization_remainder.csv"),
            row.names = FALSE)
  write.csv(remainder_table,
            file.path(out_dir, "kde_linearization_remainder.csv"),
            row.names = FALSE)
  write.csv(covariance_table, file.path(out_dir, "joint_normality_summary.csv"),
            row.names = FALSE)
  write.csv(J, file.path(out_dir, "theoretical_information.csv"))
  write.csv(Sigma, file.path(out_dir, "theoretical_asymptotic_covariance.csv"))
  write.csv(information_diagnostics,
            file.path(out_dir, "information_diagnostics.csv"),
            row.names = FALSE)
  if (nrow(bias_decomposition_table)) {
    write.csv(
      bias_decomposition_table,
      file.path(out_dir, "standardized_bias_decomposition.csv"),
      row.names = FALSE
    )
  }
  plot_marginal_qq(standardized_by_n,
                   file.path(out_dir, "marginal_standardized_QQ.pdf"))
  plot_joint_qq(whitened_by_n,
                file.path(out_dir, "joint_chisquare_QQ.pdf"))
  plot_marginal_qq(linear_standardized_by_n,
                   file.path(out_dir, "linear_benchmark_marginal_QQ.pdf"))
  plot_joint_qq(linear_whitened_by_n,
                file.path(out_dir, "linear_benchmark_joint_QQ.pdf"))
  kde_linear_standardized_plot <- lapply(
    kde_linear_standardized_by_n,
    function(z) z[complete.cases(z), , drop = FALSE]
  )
  kde_linear_whitened_plot <- lapply(
    kde_linear_whitened_by_n,
    function(z) z[complete.cases(z), , drop = FALSE]
  )
  plot_marginal_qq(kde_linear_standardized_plot,
                   file.path(out_dir, "kde_linear_marginal_QQ.pdf"))
  plot_joint_qq(kde_linear_whitened_plot,
                file.path(out_dir, "kde_linear_joint_QQ.pdf"))
  score_one_step_standardized_plot <- lapply(
    score_one_step_standardized_by_n,
    function(z) z[complete.cases(z), , drop = FALSE]
  )
  score_one_step_whitened_plot <- lapply(
    score_one_step_whitened_by_n,
    function(z) z[complete.cases(z), , drop = FALSE]
  )
  plot_marginal_qq(
    score_one_step_standardized_plot,
    file.path(out_dir, "kde_score_one_step_marginal_QQ.pdf")
  )
  plot_joint_qq(
    score_one_step_whitened_plot,
    file.path(out_dir, "kde_score_one_step_joint_QQ.pdf")
  )
  if (score_bootstrap_B > 0L) {
    bootstrap_plot <- function(z) {
      lapply(z, function(x) x[complete.cases(x), , drop = FALSE])
    }
    plot_marginal_qq(
      bootstrap_plot(score_bootstrap_full_standardized_by_n),
      file.path(out_dir, "kde_score_bootstrap_full_marginal_QQ.pdf")
    )
    plot_joint_qq(
      bootstrap_plot(score_bootstrap_full_whitened_by_n),
      file.path(out_dir, "kde_score_bootstrap_full_joint_QQ.pdf")
    )
    plot_marginal_qq(
      bootstrap_plot(score_bootstrap_shrunk_standardized_by_n),
      file.path(out_dir, "kde_score_bootstrap_shrunk_marginal_QQ.pdf")
    )
    plot_joint_qq(
      bootstrap_plot(score_bootstrap_shrunk_whitened_by_n),
      file.path(out_dir, "kde_score_bootstrap_shrunk_joint_QQ.pdf")
    )
  }
  saveRDS(list(config = list(B = B, n_values = n_values, seed = seed,
                             setting = setting,
                             beta = beta, beta_lower = beta_lower,
                             bandwidth_constant = bandwidth_constant,
                             lambda_constant = lambda_constant,
                             lambda_power = lambda_power, gamma = gamma,
                             kde_grid_n_max = kde_grid_n_max,
                             kde_points_per_bandwidth =
                               kde_points_per_bandwidth,
                             include_population_target =
                               include_population_target,
                             population_target_maxeval =
                               population_target_maxeval,
                             population_kernel_nodes =
                               population_kernel_nodes,
                             newton_hessian_step =
                               newton_hessian_step,
                             local_trust_radius =
                               local_trust_radius,
                             score_ridge_relative =
                               score_ridge_relative,
                             score_max_rootn_step =
                               score_max_rootn_step,
                             score_bootstrap_B = score_bootstrap_B,
                             score_bootstrap_min_valid_fraction =
                               score_bootstrap_min_valid_fraction),
               truth = truth, J = J, Sigma = Sigma,
               information_diagnostics = information_diagnostics,
               population_smoothed_target = population_target_table,
               standardized_bias_decomposition =
                 bias_decomposition_table,
               standard_deviation_comparison =
                 standard_deviation_comparison_table,
               marginal_sd_coverage = marginal_sd_coverage_table,
               raw_error = raw_by_n, standardized = standardized_by_n,
               whitened = whitened_by_n, fit_diagnostics = diagnostics_table,
               failure_diagnostics = failure_table,
               linear_raw_error = linear_raw_by_n,
               linear_standardized = linear_standardized_by_n,
               linear_whitened = linear_whitened_by_n,
               kde_linear_raw_error = kde_linear_raw_by_n,
               kde_linear_standardized = kde_linear_standardized_by_n,
               kde_linear_whitened = kde_linear_whitened_by_n,
               fisher_newton_raw_error = fisher_newton_raw_by_n,
               fisher_newton_standardized =
                 fisher_newton_standardized_by_n,
               observed_newton_raw_error = observed_newton_raw_by_n,
               observed_newton_standardized =
                 observed_newton_standardized_by_n,
               kde_score_one_step_raw_error = score_one_step_raw_by_n,
               kde_score_one_step_standardized =
                 score_one_step_standardized_by_n,
               kde_score_one_step_whitened = score_one_step_whitened_by_n,
               kde_score_one_step_marginal_summary =
                 score_one_step_summary_table,
               kde_score_one_step_joint_summary =
                 score_one_step_joint_table,
               kde_score_one_step_diagnostics =
                 score_one_step_diagnostics_table,
               kde_score_one_step_remainder =
                 score_one_step_remainder_table,
               kde_score_bootstrap_full_raw_error =
                 score_bootstrap_full_raw_by_n,
               kde_score_bootstrap_full_standardized =
                 score_bootstrap_full_standardized_by_n,
               kde_score_bootstrap_full_whitened =
                 score_bootstrap_full_whitened_by_n,
               kde_score_bootstrap_full_marginal_summary =
                 score_bootstrap_full_summary_table,
               kde_score_bootstrap_full_joint_summary =
                 score_bootstrap_full_joint_table,
               kde_score_bootstrap_shrunk_raw_error =
                 score_bootstrap_shrunk_raw_by_n,
               kde_score_bootstrap_shrunk_standardized =
                 score_bootstrap_shrunk_standardized_by_n,
               kde_score_bootstrap_shrunk_whitened =
                 score_bootstrap_shrunk_whitened_by_n,
               kde_score_bootstrap_shrunk_marginal_summary =
                 score_bootstrap_shrunk_summary_table,
               kde_score_bootstrap_shrunk_joint_summary =
                 score_bootstrap_shrunk_joint_table,
               kde_score_bootstrap_diagnostics =
                 score_bootstrap_diagnostics_table,
               fisher_newton_marginal_summary =
                 fisher_newton_summary_table,
               observed_newton_marginal_summary =
                 observed_newton_summary_table,
               newton_decomposition_summary = newton_comparison_table,
               newton_target_diagnostics = newton_diagnostics_table,
               linear_marginal_summary = linear_summary_table,
               linear_joint_summary = linear_joint_table,
               kde_linear_marginal_summary = kde_linear_summary_table,
               kde_linear_joint_summary = kde_linear_joint_table,
               linearization_remainder = remainder_table,
               marginal_summary = summary_table,
               joint_summary = covariance_table, fits = all_results),
          file.path(out_dir, "scenario1c_normality_results.rds"))

  cat("\nCompleted. Results written to:\n", normalizePath(out_dir), "\n")
  print(diagnostics_table, row.names = FALSE)
  print(covariance_table, row.names = FALSE)
  invisible(list(marginal = summary_table, joint = covariance_table,
                 diagnostics = diagnostics_table,
                 population_target = population_target_table,
                 bias_decomposition = bias_decomposition_table,
                 standard_deviation = standard_deviation_comparison_table,
                 newton_decomposition = newton_comparison_table,
                 newton_diagnostics = newton_diagnostics_table,
                 J = J, Sigma = Sigma))
}

parse_cli <- function(args) {
  cfg <- list(B = 200L, n_values = c(1000L, 2000L, 5000L),
              ncores = 6L, quick = FALSE, out_dir = NULL,
              setting = "s1c", seed = NULL,
              beta = NULL, bandwidth_constant = 0.5,
              kde_grid_n_max = 131072L,
              kde_points_per_bandwidth = 10,
              include_population_target = TRUE,
              population_target_maxeval = 2000L,
              population_kernel_nodes = 401L,
              newton_hessian_step = 1e-4,
              local_trust_radius = 8,
              score_ridge_relative = 1e-3,
              score_max_rootn_step = 6,
              score_bootstrap_B = 0L,
              score_bootstrap_min_valid_fraction = 0.7)
  for (a in args) {
    if (a == "--quick") cfg$quick <- TRUE
    else if (a == "--sequential") cfg$ncores <- 1L
    else if (a == "--skip-population-target") {
      cfg$include_population_target <- FALSE
    }
    else if (grepl("^--B=", a)) cfg$B <- as.integer(sub("^--B=", "", a))
    else if (grepl("^--setting=", a)) {
      cfg$setting <- tolower(sub("^--setting=", "", a))
    } else if (grepl("^--seed=", a)) {
      cfg$seed <- as.integer(sub("^--seed=", "", a))
    }
    else if (grepl("^--n=", a)) {
      cfg$n_values <- as.integer(strsplit(sub("^--n=", "", a), ",", fixed = TRUE)[[1]])
    } else if (grepl("^--cores=", a)) {
      cfg$ncores <- as.integer(sub("^--cores=", "", a))
    } else if (grepl("^--beta=", a)) {
      cfg$beta <- as.numeric(sub("^--beta=", "", a))
    } else if (grepl("^--bw-constant=", a)) {
      cfg$bandwidth_constant <-
        as.numeric(sub("^--bw-constant=", "", a))
    } else if (grepl("^--grid-max=", a)) {
      cfg$kde_grid_n_max <- as.integer(sub("^--grid-max=", "", a))
    } else if (grepl("^--points-per-h=", a)) {
      cfg$kde_points_per_bandwidth <-
        as.numeric(sub("^--points-per-h=", "", a))
    } else if (grepl("^--population-maxeval=", a)) {
      cfg$population_target_maxeval <-
        as.integer(sub("^--population-maxeval=", "", a))
    } else if (grepl("^--population-kernel-nodes=", a)) {
      cfg$population_kernel_nodes <-
        as.integer(sub("^--population-kernel-nodes=", "", a))
    } else if (grepl("^--newton-hessian-step=", a)) {
      cfg$newton_hessian_step <-
        as.numeric(sub("^--newton-hessian-step=", "", a))
    } else if (grepl("^--local-radius=", a)) {
      cfg$local_trust_radius <-
        as.numeric(sub("^--local-radius=", "", a))
    } else if (grepl("^--score-ridge=", a)) {
      cfg$score_ridge_relative <-
        as.numeric(sub("^--score-ridge=", "", a))
    } else if (grepl("^--score-max-step=", a)) {
      cfg$score_max_rootn_step <-
        as.numeric(sub("^--score-max-step=", "", a))
    } else if (grepl("^--score-bootstrap=", a)) {
      cfg$score_bootstrap_B <-
        as.integer(sub("^--score-bootstrap=", "", a))
    } else if (grepl("^--score-bootstrap-min-valid=", a)) {
      cfg$score_bootstrap_min_valid_fraction <-
        as.numeric(sub("^--score-bootstrap-min-valid=", "", a))
    } else if (grepl("^--out=", a)) {
      cfg$out_dir <- sub("^--out=", "", a)
    } else stop("Unknown command-line option: ", a)
  }
  if (cfg$quick) {
    cfg$B <- 2L
    cfg$n_values <- 500L
    cfg$ncores <- 1L
  }
  cfg
}

if (sys.nframe() == 0L) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  call_args <- list(
    B = cli$B, n_values = cli$n_values, ncores = cli$ncores,
    setting = cli$setting, seed = cli$seed,
    beta = cli$beta, bandwidth_constant = cli$bandwidth_constant,
    kde_grid_n_max = cli$kde_grid_n_max,
    kde_points_per_bandwidth = cli$kde_points_per_bandwidth,
    include_population_target = cli$include_population_target,
    population_target_maxeval = cli$population_target_maxeval,
    population_kernel_nodes = cli$population_kernel_nodes,
    newton_hessian_step = cli$newton_hessian_step,
    local_trust_radius = cli$local_trust_radius,
    score_ridge_relative = cli$score_ridge_relative,
    score_max_rootn_step = cli$score_max_rootn_step,
    score_bootstrap_B = cli$score_bootstrap_B,
    score_bootstrap_min_valid_fraction =
      cli$score_bootstrap_min_valid_fraction
  )
  if (!is.null(cli$out_dir)) call_args$out_dir <- cli$out_dir
  if (isTRUE(cli$quick)) {
    call_args$nstart <- 1L
    call_args$max_iter <- 20L
    call_args$maxeval <- 80L
    call_args$information_grid_n <- 20001L
    call_args$population_target_maxeval <- 100L
  }
  do.call(run_scenario1c_normality, call_args)
}
